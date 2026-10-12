within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setSat_T_Newton
  "saturated liquid and vapour density at temperature T: Newton's method of Akasaka (2008) with backtracking"
  input Temperature T;
  input Real delta_liq_start(unit="1") "start value for the reduced liquid density";
  input Real delta_vap_start(unit="1") "start value for the reduced vapour density";
  output Boolean converged=false;
  output Real delta_liq(unit="1") "reduced saturated liquid density";
  output Real delta_vap(unit="1") "reduced saturated vapour density";
  output Real RES[2] "residuals of the dimensionless pressure J and Gibbs energy K";

protected
  constant Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  constant Density d_reducing=helmholtzCoefficients.d_reducing "reducing density of the EoS";
  constant Real delta_crit=d_crit/d_reducing "critical density, reduced with the EoS reducing density";
  constant Real delta_liq_max=1.1*fluidLimits.DMAX/d_reducing "upper bound for reduced liquid density";
  constant Real delta_vap_min=fluidLimits.DMIN/d_reducing "lower bound for reduced vapour density";
  constant Real tolerance=1e-9 "tolerance for the residuals of J (relative) and K (absolute)";
  constant Real tolerance_step=1e-7 "tolerance for the Newton step, relative to the density gap delta_liq-delta_vap";
  constant Integer iter_max=100;
  constant Integer backtrack_max=30;
  Real tau(unit="1")=helmholtzCoefficients.T_reducing/T "inverse reduced temperature";

  Real r_liq; Real rd_liq; Real rdd_liq;
  Real r_vap; Real rd_vap; Real rdd_vap;
  Real J_liq; Real J_vap; Real K_liq; Real K_vap;
  Real J_liq_delta; Real J_vap_delta; Real K_liq_delta; Real K_vap_delta;
  Real Jacobian[2,2] "Jacobian matrix";
  Real NS[2] "Newton step vector";
  Real lambda "step length, halved until the step is acceptable";
  Real delta_liq_new; Real delta_vap_new;
  Real RES_new[2];
  Boolean accepted;
  Integer iter=0;
  Integer k;

algorithm
  delta_liq := min(max(delta_liq_start, delta_crit*(1+1e-12)), delta_liq_max);
  delta_vap := min(max(delta_vap_start, delta_vap_min), delta_crit*(1-1e-12));

  r_liq := EoS.f_r(tau=tau, delta=delta_liq);
  rd_liq := EoS.f_rd(tau=tau, delta=delta_liq);
  rdd_liq := EoS.f_rdd(tau=tau, delta=delta_liq);
  r_vap := EoS.f_r(tau=tau, delta=delta_vap);
  rd_vap := EoS.f_rd(tau=tau, delta=delta_vap);
  rdd_vap := EoS.f_rdd(tau=tau, delta=delta_vap);
  J_liq := delta_liq*(1 + delta_liq*rd_liq);
  J_vap := delta_vap*(1 + delta_vap*rd_vap);
  K_liq := delta_liq*rd_liq + r_liq + log(delta_liq);
  K_vap := delta_vap*rd_vap + r_vap + log(delta_vap);
  J_liq_delta := 1 + 2*delta_liq*rd_liq + delta_liq*delta_liq*rdd_liq;
  J_vap_delta := 1 + 2*delta_vap*rd_vap + delta_vap*delta_vap*rdd_vap;
  K_liq_delta := 2*rd_liq + delta_liq*rdd_liq + 1/delta_liq;
  K_vap_delta := 2*rd_vap + delta_vap*rdd_vap + 1/delta_vap;
  RES := {J_vap-J_liq, K_vap-K_liq};

  accepted := true;
  while (iter<iter_max) and accepted and not converged loop
    iter := iter+1;
    Jacobian := [-J_liq_delta, J_vap_delta;
                 -K_liq_delta, K_vap_delta];
    NS := -Modelica.Math.Matrices.solve(Jacobian, RES);

    // backtracking: halve the step until both phases stay on their side of the critical density,
    // remain mechanically stable (dJ/ddelta>0) and the residual does not increase
    lambda := 1;
    accepted := false;
    k := 0;
    while (k<backtrack_max) and not accepted loop
      k := k+1;
      delta_liq_new := delta_liq + lambda*NS[1];
      delta_vap_new := delta_vap + lambda*NS[2];
      if (delta_liq_new>delta_crit) and (delta_liq_new<delta_liq_max) and (delta_vap_new>delta_vap_min) and (delta_vap_new<delta_crit) then
        r_liq := EoS.f_r(tau=tau, delta=delta_liq_new);
        rd_liq := EoS.f_rd(tau=tau, delta=delta_liq_new);
        rdd_liq := EoS.f_rdd(tau=tau, delta=delta_liq_new);
        r_vap := EoS.f_r(tau=tau, delta=delta_vap_new);
        rd_vap := EoS.f_rd(tau=tau, delta=delta_vap_new);
        rdd_vap := EoS.f_rdd(tau=tau, delta=delta_vap_new);
        J_liq_delta := 1 + 2*delta_liq_new*rd_liq + delta_liq_new*delta_liq_new*rdd_liq;
        J_vap_delta := 1 + 2*delta_vap_new*rd_vap + delta_vap_new*delta_vap_new*rdd_vap;
        J_liq := delta_liq_new*(1 + delta_liq_new*rd_liq);
        J_vap := delta_vap_new*(1 + delta_vap_new*rd_vap);
        K_liq := delta_liq_new*rd_liq + r_liq + log(delta_liq_new);
        K_vap := delta_vap_new*rd_vap + r_vap + log(delta_vap_new);
        RES_new := {J_vap-J_liq, K_vap-K_liq};
        accepted := (J_liq_delta>0) and (J_vap_delta>0) and (RES_new*RES_new <= RES*RES);
      end if;
      if not accepted then
        lambda := lambda/2;
      end if;
    end while;

    if not accepted then
      // no step reduces the residual any further; close to the critical point this happens when the residuals
      // are at roundoff level, before the step criterion is met: this is the best attainable solution
      converged := abs(RES[1])<=1e-13*J_vap and abs(RES[2])<=1e-13;
    end if;

    if accepted then
      delta_liq := delta_liq_new;
      delta_vap := delta_vap_new;
      RES := RES_new;
      K_liq_delta := 2*rd_liq + delta_liq*rdd_liq + 1/delta_liq;
      K_vap_delta := 2*rd_vap + delta_vap*rdd_vap + 1/delta_vap;
      // residuals: J relative to J_vap, or for the (almost incompressible) liquid at low pressure measured as density error;
      // close to the critical point small residuals do not imply an accurate solution, therefore the Newton step
      // must also be small compared to the density gap
      converged := (abs(RES[1])<=tolerance*J_vap or abs(RES[1]/J_liq_delta)<=tolerance*delta_liq)
               and abs(RES[2])<=tolerance
               and max(abs(lambda*NS[1]), abs(lambda*NS[2]))<=tolerance_step*(delta_liq-delta_vap);
    end if;
  end while;
end setSat_T_Newton;
