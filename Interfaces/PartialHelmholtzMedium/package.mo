within HelmholtzMedia.Interfaces;
partial package PartialHelmholtzMedium 
  extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(
    onePhase=false,
    singleState=false,
    smoothModel=true,
    DipoleMoment(min=0, max=5),
    AbsolutePressure(min=Modelica.Constants.small, max=1e12),
    SpecificEntropy(min=-Modelica.Constants.inf, max=Modelica.Constants.inf));

import HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Types.*;

  constant FluidLimits fluidLimits;

  constant EoS.HelmholtzCoefficients helmholtzCoefficients;

  constant Transport.ThermalConductivityCoefficients thermalConductivityCoefficients;

  constant Transport.DynamicViscosityCoefficients dynamicViscosityCoefficients;

  constant Transport.SurfaceTensionCoefficients surfaceTensionCoefficients;

  constant Ancillary.AncillaryCoefficients ancillaryCoefficients;

  // constant IndependentVariables independentVariables=IndependentVariables.dTX;


  redeclare record extends ThermodynamicState(phase(start=0))
    // inherits phase integer
    Temperature T "Temperature of medium";
    AbsolutePressure p "Absolute pressure of medium";
    Density d "Density of medium";
    SpecificEnergy u "Specific inner energy of medium";
    SpecificEnthalpy h "Specific enthalpy of medium";
    SpecificEntropy s "Specific entropy of medium";
  end ThermodynamicState;


  redeclare record extends SaturationProperties
    // inherits Tsat and psat
    ThermodynamicState liq;
    ThermodynamicState vap;
  end SaturationProperties;


  redeclare model extends BaseProperties
  "Base properties (p, d, T, h, u, s) of a medium"

    SpecificEntropy s;

  equation
    MM = fluidConstants[1].molarMass;
    R = Modelica.Constants.R/MM;

  // use functions to calculate properties
    d = density_ph(p=p, h=h);
    T = temperature_ph(p=p, h=h);
    // s = specificEntropy_ph(p=p, h=h);

    // d = density_pT(p=p, T=T);
    // p =  pressure_dT(d=d, T=T);
    // h =  specificEnthalpy_dT(d=d, T=T);
    s =  specificEntropy_dT(d=d, T=T);

  // calculate u
    u = h - p/d;

  // set SaturationProperties
    sat = setSat_p(p=p);

  // connect state with BaseProperties
    state.d = d;
    state.T = T;
    state.p = p;
    state.h = h;
    state.s = s;
    state.u = u;
    state.phase =if ((p < fluidConstants[1].criticalPressure) and (p > fluidConstants[1].triplePointPressure) and (h > sat.liq.h) and (h < sat.vap.h)) then 2 else 1;

  end BaseProperties;


  redeclare function setSat_T
  "iterative calculation of saturation properties from EoS with Newton-Raphson algorithm"
    // does not extend, because base class already defines an algorithm
    input Temperature T;
    output SaturationProperties sat;

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real tau(unit="1") "inverse reduced temperature";

    Real delta_liq(unit="1", min=0);
    Real delta_vap(unit="1", min=0);
    Real J_liq;
    Real J_liq_delta;
    Real J_vap;
    Real J_vap_delta;
    Real Delta_J;
    Real K_liq;
    Real K_liq_delta;
    Real K_vap;
    Real K_vap_delta;
    Real Delta_K;
    Real det "determinant of Jacobi matrix";
    Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
    Integer iter=0;
    Real tolerance=1e-9 "Tolerance for sum of Delta_J and Delta_K";

  algorithm
    // Modelica.Utilities.Streams.print("setSat_T: T="+String(T),"printlog.txt");

  if ((T>=T_trip) and (T<T_crit)) then
    tau := T_crit/T;

    // calculate guess values for reduced density delta
    delta_liq := Ancillary.bubbleDensity_T(T=T)/d_crit;
    delta_vap := Ancillary.dewDensity_T(T=T)/d_crit;

    // pressure difference liquid-vapor
    J_liq := delta_liq*(1 + delta_liq*EoS.f_rd(tau=tau, delta=delta_liq));
    J_vap := delta_vap*(1 + delta_vap*EoS.f_rd(tau=tau, delta=delta_vap));
    Delta_J := (J_vap-J_liq);

    // Gibbs energy difference liquid-vapor
    K_liq := delta_liq*EoS.f_rd(tau=tau, delta=delta_liq) + EoS.f_r(tau=tau, delta=delta_liq) + log(delta_liq);
    K_vap := delta_vap*EoS.f_rd(
                            tau=tau, delta=delta_vap) + EoS.f_r(tau=tau, delta=delta_vap) + log(delta_vap);
    Delta_K := (K_vap-K_liq);

    while (abs(Delta_J) + abs(Delta_K) > tolerance) loop
      // Modelica.Utilities.Streams.print(" ", "printlog.txt");
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("delta_liq=" + String(delta_liq) + " and delta_vap=" + String(delta_vap), "printlog.txt");
      // Modelica.Utilities.Streams.print("Delta_J=" + String(Delta_J) + " and Delta_K=" + String(Delta_K), "printlog.txt");
      iter := iter+1;

      // calculate gradients
      J_liq_delta := 1 + 2*delta_liq*EoS.f_rd(tau=tau, delta=delta_liq) + delta_liq^2*EoS.f_rdd(tau=tau, delta=delta_liq);
      J_vap_delta := 1 + 2*delta_vap*EoS.f_rd(tau=tau, delta=delta_vap) + delta_vap^2*EoS.f_rdd(tau=tau, delta=delta_vap);
      K_liq_delta := 2*EoS.f_rd(tau=tau, delta=delta_liq) + delta_liq*EoS.f_rdd(tau=tau, delta=delta_liq) + 1/delta_liq;
      K_vap_delta := 2*EoS.f_rd(tau=tau, delta=delta_vap) + delta_vap*EoS.f_rdd(tau=tau, delta=delta_vap) + 1/delta_vap;

      // calculate determinant of Jacobi matrix
      det := J_vap_delta*K_liq_delta - J_liq_delta*K_vap_delta;

      // calculate better values for reduced density delta
      delta_liq := delta_liq + gamma/det*((K_vap - K_liq)*J_vap_delta - (J_vap - J_liq)*K_vap_delta);
      delta_vap := delta_vap + gamma/det*((K_vap - K_liq)*J_liq_delta - (J_vap - J_liq)*K_liq_delta);

      // check bounds
      delta_liq := max(Modelica.Constants.small, delta_liq);
      delta_vap := max(Modelica.Constants.small, delta_vap);

      // calculate new Delta_J and Delta_K
      J_liq := delta_liq*(1 + delta_liq*EoS.f_rd(tau=tau, delta=delta_liq));
      J_vap := delta_vap*(1 + delta_vap*EoS.f_rd(tau=tau, delta=delta_vap));
      Delta_J := (J_vap-J_liq);

      K_liq := delta_liq*EoS.f_rd(tau=tau, delta=delta_liq) + EoS.f_r(tau=tau, delta=delta_liq) + log(delta_liq);
      K_vap := delta_vap*EoS.f_rd(tau=tau, delta=delta_vap) + EoS.f_r(tau=tau, delta=delta_vap) + log(delta_vap);
      Delta_K := (K_vap-K_liq);
    end while;
    // Modelica.Utilities.Streams.print("setSat_T total iteration steps " + String(iter), "printlog.txt");

    sat.Tsat := T;
    sat.liq := setState_dTX(d=delta_liq*d_crit, T=T, phase=1);
    sat.vap := setState_dTX(d=delta_vap*d_crit, T=T, phase=1);
    sat.psat := sat.liq.p;

  elseif (T>=T_crit) then
    // assert(T <= T_crit, "setSat_T error: Temperature is higher than critical temperature");
    // above critical temperature, no stable two-phase state exists
    // anyway, it is possible to extend the vapour-pressure curve into this region
    // this can happen when called from BaseProperties
    // one possibility is use the state where ds/dT=max or ds/dp=max or dcp/dT=max or dcp/dp=max
    // here a very simple approximation is used by just setting d=d_crit
    sat.Tsat := T;
    sat.liq := setState_dTX(d=d_crit, T=T, phase=1);
    sat.vap := setState_dTX(d=d_crit, T=T, phase=1);
    sat.psat := sat.liq.p;
  else
    // assert(T >= T_trip, "setSat_T error: Temperature is lower than triple-point temperature");
    // T<T_trip: this does not make sense: if T is below the triple temperature, the medium is solid, not fluid
    // anyway, during initialization (at time=0) T=0 may happen
    // density values are extrapolated linearly, fantasy values are returned
    sat.Tsat := max(T, Modelica.Constants.small);
    delta_liq := (T_trip/sat.Tsat)*Ancillary.bubbleDensity_T(T=T_trip)/d_crit;
    delta_vap := (sat.Tsat/T_trip)*Ancillary.dewDensity_T(T=T_trip)/d_crit;

    sat.liq := setState_dTX(d=delta_liq*d_crit, T=T, phase=1);
    sat.vap := setState_dTX(d=delta_vap*d_crit, T=T, phase=1);
    sat.psat := sat.liq.p;
  end if;

    annotation (Documentation(info="<html>
This function iteratively determines the saturation state  for a given temperature 
by varying the density of saturated liquid and saturated vapor 
with a Newton-Raphson approach for simultaneous equations.

<dl>
<dt> Ryo Akasaka:</dt>
<dd> <b>A reliable and useful method to determine the saturation state from Helmholtz energy equations of state</b>.<br>
     Journal of Thermal Science and Technology 3 (3) , 442-451.<br>
     DOI: <a href=\"http://dx.doi.org/10.1299/jtst.3.442\">10.1299/jtst.3.442</a>
</dd>
</dl>
</html>"));
  end setSat_T;


  redeclare function setSat_p
  "iterative calculation of saturation properties from EoS for a given pressure"
    input AbsolutePressure p;
    output SaturationProperties sat;

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    constant Density dv_trip = Ancillary.dewDensity_T(T_trip);
    constant Density dl_trip = Ancillary.bubbleDensity_T(T_trip);

    EoS.HelmholtzDerivs fl;
    EoS.HelmholtzDerivs fv;

    Real RES[3] "residual function vector";
    Real RSS "residual sum of squares";
    Real Jacobian[3,3] "Jacobian matrix";
    Real NS[3] "Newton step vector";

    constant Real lambda(min=0.1,max=1) = 1 "convergence speed, default=1";
    constant Real tolerance=1e-5 "tolerance for RSS";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("setSat_p: p=" + String(p), "printlog.txt");

  if ((p>p_trip) and (p<p_crit)) then
    // calculate start values
    // sat.Tsat  := 1/(1/T_crit - (1/T_trip-1/T_crit)/log(p_crit/p_trip)*log(p/p_crit));
    // at lower p the difference between dl and dv is bigger
    sat.Tsat := 0.99*Ancillary.saturationTemperature_p(p=p);
    sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
    sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);

    // calculate residuals: RES=calc-input
    fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
    fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
    RES := {EoS.p(fl)-p, EoS.p(fv)-p, EoS.g(fl)-EoS.g(fv)};
    RSS := RES*RES/2;

    while (RSS>tolerance) and (iter<iter_max) loop
      iter := iter+1;

      // calculate Jacobian matrix and Newton Step vector
      Jacobian := [+EoS.dpdT(fl), +0,            +EoS.dpTd(fl);
                   +0,            +EoS.dpdT(fv), +EoS.dpTd(fv);
                   +EoS.dgdT(fl), -EoS.dgdT(fv), EoS.dgTd(fl)-EoS.dgTd(fv)];
      NS := -Modelica.Math.Matrices.solve(Jacobian,RES);

      // calculate better sat.liq.d, sat.vap.d and sat.Tsat
      sat.liq.d := sat.liq.d + lambda*NS[1];
      sat.vap.d := sat.vap.d + lambda*NS[2];
      sat.Tsat  := sat.Tsat  + lambda*NS[3];

      // check bounds
      sat.liq.d := max(sat.liq.d, 0.98*d_crit);
      sat.liq.d := min(sat.liq.d, 1.02*dl_trip);
      sat.vap.d := max(sat.vap.d, 0.98*dv_trip);
      sat.vap.d := min(sat.vap.d, 1.02*d_crit);
      sat.Tsat  := max(sat.Tsat,  0.98*T_trip);
      sat.Tsat  := min(sat.Tsat,  1.02*T_crit);

      // calculate new residuals: RES=calc-input
      fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
      fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
      RES := {EoS.p(fl)-p, EoS.p(fv)-p, EoS.g(fl)-EoS.g(fv)};
      RSS := RES*RES/2;
    end while;
    // if verbose then Modelica.Utilities.Streams.print("setSat_p total iteration steps " + String(iter), "printlog.txt"); end if;
    // Modelica.Utilities.Streams.print("setSat_p total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setSat_p did not converge, input was p=" + String(p) + "; the remaining residuals are RES_pl=" + String(RES[1]) + " and RES_pv=" + String(RES[2]) + " and RES_g=" + String(RES[3]));
    // check bounds, more strict
    sat.liq.d := max(sat.liq.d, d_crit);
    sat.liq.d := min(sat.liq.d, dl_trip);
    sat.vap.d := max(sat.vap.d, dv_trip);
    sat.vap.d := min(sat.vap.d, d_crit);
    sat.Tsat  := max(sat.Tsat,  T_trip);
    sat.Tsat  := min(sat.Tsat,  T_crit);

  elseif (p>=p_crit) then
    // assert(p <= p_crit, "setSat_p error: pressure is higher than critical pressure");
    // above critical pressure, no stable two-phase state exists
    // anyway, it is possible to extend the vapour-pressure curve into this region
    // this can happen when called from BaseProperties
    // one possibility is to use the state where ds/dT=max or ds/dp=max or dcp/dT=max or dcp/dp=max
    // here, critical values are returned
    sat.psat  := p_crit;
    sat.Tsat  := T_crit;
    sat.liq.d := d_crit;
    sat.vap.d := d_crit;
  elseif (p<=p_trip) then
    sat.psat  := p;
    sat.Tsat  := T_trip;
    sat.liq.d := dl_trip;
    sat.vap.d := dv_trip;
  else
    assert(false, "setSat_p: this should not happen, check p");
  end if;

    sat.psat := p;
    sat.liq := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);

  end setSat_p;


  redeclare function extends setBubbleState
  "returns bubble ThermodynamicState from given saturation properties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat, input phase and output state
  algorithm
    state := sat.liq;
  end setBubbleState;


  redeclare function extends setDewState
  "returns dew ThermodynamicState from given saturation properties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat, input phase and output state
  algorithm
    state := sat.vap;
  end setDewState;


  redeclare function extends setState_dTX
  "Return thermodynamic state as function of (d, T)"

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;

    EoS.HelmholtzDerivs f;
    SaturationProperties sat;
    MassFraction x "vapour quality";

  algorithm
    state.phase := phase;

    if (state.phase == 0) then
      //phase unknown, check phase first
      if ((T>=T_trip) and (T<T_crit)) then
        // two-phase possible, do simple density check
        // Modelica.Utilities.Streams.print("setState_dT: dliq=" + String(bubbleDensity_T_ANC(T=T)) + " dvap=" + String(dewDensity_T_ANC(T=T)) + ", simple check only");
        if ((d > 1.05*Ancillary.bubbleDensity_T(T=T)) or (d < 0.98*Ancillary.dewDensity_T(T=T))) then
          state.phase := 1;
        else
          // Modelica.Utilities.Streams.print("setState_dT: d=" + String(d) + " T=" + String(T) + ", two-phase state or close to it");
          // get saturation properties from EoS, use Tsat as starting value
          sat := setSat_T(T=T);
          if ((d < sat.liq.d) and (d > sat.vap.d)) then
            state.phase := 2;
          else
            state.phase := 1;
          end if;
        end if;
      else
        // T>=T_crit
        state.phase := 1;
      end if;
    elseif (state.phase == 2) then
      assert(T <= T_crit, "setState_dTX_error: Temperature is higher than critical temperature");
      sat := setSat_T(T=T);
      assert(d >= sat.vap.d, "setState_dTX_error: density is lower than saturated vapor density: this is single phase vapor");
      assert(d <= sat.liq.d, "setState_dTX_error: density is higher than saturated liquid density: this is single phase liquid");
    end if;

    state.d := d;
    state.T := T;
    if (state.phase == 2) then
      // force two-phase
      x := (1.0/d - 1.0/sat.liq.d)/(1.0/sat.vap.d - 1.0/sat.liq.d);
      state.p := sat.psat;
      state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
    else
      // force single-phase
      f := EoS.setHelmholtzDerivsFirst(d=d,T=T);
      state.p := EoS.p(f=f);
      state.s := EoS.s(f=f);
      state.h := EoS.h(f=f);
      state.u := EoS.u(f=f);
    end if;

  end setState_dTX;


  redeclare function setState_Tx
  "Return thermodynamic state as function of (T, x)"
    input Temperature T "Temperature";
    input MassFraction x "Vapour quality";
    output ThermodynamicState state "Thermodynamic state record";

protected
    SaturationProperties sat=setSat_T(T=T);

  algorithm
    state.phase := 2;
    state.p := sat.psat;
    state.T := sat.Tsat;
    state.d := 1.0/(1.0/sat.liq.d + x*(1.0/sat.vap.d - 1.0/sat.liq.d));
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
    state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
  end setState_Tx;


  redeclare function setState_px
  "Return thermodynamic state as function of (p, x)"
    input AbsolutePressure p "Pressure";
    input MassFraction x "Vapour quality";
    output ThermodynamicState state "Thermodynamic state record";

protected
    SaturationProperties sat=setSat_p(p=p);

  algorithm
    state.phase := 2;
    state.p := sat.psat;
    state.T := sat.Tsat;
    state.d := 1.0/(1.0/sat.liq.d + x*(1.0/sat.vap.d - 1.0/sat.liq.d));
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
    state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
  end setState_px;


  redeclare function extends setState_pTX
  "Return thermodynamic state as function of (p, T)"

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    EoS.HelmholtzDerivs f(T=T);
    SaturationProperties sat;

    Density d_min;
    Density d_max;
    Density d_iter;
    AbsolutePressure RES_p;
    AbsolutePressure RES_min;
    AbsolutePressure RES_max;
    DerPressureByDensity dpdT "(dp/dd)@T=const";
    Real gamma(min=0.1,max=1) = 1 "convergence speed, default=1";
    constant Real tolerance=1e-9 "relativ tolerance for RES_p";
    Integer iter = 0;
    constant Integer iter_max=200;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("p=" + String(p) + " and T=" + String(T), "printlog.txt");
    assert(phase <> 2, "setState_pTX_error: pressure and temperature are not independent variables in two-phase state");
    state.phase := 1;

    if (T < T_crit) then
      // determine p_sat
      sat.psat := Ancillary.saturationPressure_T(T=T);
      if (p > 1.02*sat.psat) then
        sat.liq.d := Ancillary.bubbleDensity_T(T=T);
      elseif (p < 0.98*sat.psat) then
        sat.vap.d := Ancillary.dewDensity_T(T=T);
      else
        // Modelica.Utilities.Streams.print("close to saturation boundary, get saturation properties from EoS", "printlog.txt");
        sat := setSat_T(T=T);
      end if;
      // Modelica.Utilities.Streams.print("sat.psat=" + String(sat.psat), "printlog.txt");

      // determine region
      if (p > sat.psat) then
        // Modelica.Utilities.Streams.print("single-phase liquid region", "printlog.txt");
        d_min  := 0.98*sat.liq.d;
        d_max  := 1.1*fluidLimits.DMAX;
        // d_iter := 1.02*sat.liq.d;
        d_iter := Ancillary.density_pT_Soave(T=T, p=p, psat=sat.psat);
      elseif (p < sat.psat) then
        // Modelica.Utilities.Streams.print("single-phase vapor region", "printlog.txt");
        d_min  := fluidLimits.DMIN;
        d_max  := 1.02*sat.vap.d;
        // d_iter := p/(R*T);
        d_iter := Ancillary.density_pT_Soave(T=T, p=p, psat=sat.psat);
      else
        // Modelica.Utilities.Streams.print("two-phase region", "printlog.txt");
        assert(p <> sat.psat, "setState_pTX_error: pressure equals saturation pressure");
      end if;
    else
      // Modelica.Utilities.Streams.print("single-phase super-critical region", "printlog.txt");
      d_min  := fluidLimits.DMIN;
      d_max  := 1.1*fluidLimits.DMAX;
      // d_iter := p/(R*T);
      d_iter := Ancillary.density_pT_Soave(T=T, p=p, psat=sat.psat);
    end if;

    // Modelica.Utilities.Streams.print("phase and region determined, d_min=" + String(d_min) + ", d_max=" + String(d_max) + " and d_iter=" + String(d_iter), "printlog.txt");

    // check bounds, ideal gas and Soave are not very accurate
    d_iter := max(d_iter, d_min);
    d_iter := min(d_iter, d_max);

    // Modelica.Utilities.Streams.print("initialize bisection", "printlog.txt");
    // min
    f.d := d_min;
    f.delta := f.d/d_crit;
    f.rd  := EoS.f_rd(delta=f.delta, tau=f.tau);
    RES_min := EoS.p(f) - p;
    // max
    f.d := d_max;
    f.delta := f.d/d_crit;
    f.rd  := EoS.f_rd(delta=f.delta, tau=f.tau);
    RES_max := EoS.p(f) - p;
    // iter
    f.d := d_iter;
    f.delta := f.d/d_crit;
    f.rd  := EoS.f_rd(delta=f.delta, tau=f.tau);
    RES_p := EoS.p(f) - p;

    assert((RES_min*RES_max<0), "setState_pTX: d_min and d_max did not bracket the root");
    // thighten the bounds
    // opposite sign brackets the root
    if (RES_p*RES_min<0) then
      d_max := d_iter;
      RES_max := RES_p;
    elseif (RES_p*RES_max<0) then
      d_min := d_iter;
      RES_min := RES_p;
    end if;

    while ((abs(RES_p/p) > tolerance) and (iter<iter_max)) loop
      iter := iter+1;
      // gamma := iter/(iter+1);

      // calculate gradient with respect to density
      f.rdd := EoS.f_rdd(delta=f.delta, tau=f.tau);
      dpdT := EoS.dpdT(f);

      // print for Newton debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter) + ", current d_iter=" + String(d_iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("RES_p=" + String(RES_p) + " and dpdT=" + String(dpdT), "printlog.txt");

      // calculate better d_iter using Newton
      d_iter := d_iter - gamma/dpdT*RES_p;

      // check bounds, if out of bounds use bisection
      if (d_iter<d_min) or (d_iter>d_max) then
        // Modelica.Utilities.Streams.print("d_iter out of bounds, fallback to bisection method, step=" + String(iter) + ", d_iter=" + String(d_iter), "printlog.txt");
        d_iter := (d_min+d_max)/2;
      end if;

      // set necessary parts of f, then calculate new RES_p
      f.d := d_iter;
      f.delta := f.d/d_crit;
      f.rd  := EoS.f_rd(delta=f.delta, tau=f.tau);
      RES_p := EoS.p(f) - p;

      // thighten the bounds
      // opposite sign brackets the root
      if (RES_p*RES_min<0) then
        d_max := d_iter;
        RES_max := RES_p;
      elseif (RES_p*RES_max<0) then
        d_min := d_iter;
        RES_min := RES_p;
      end if;
      // Modelica.Utilities.Streams.print("d_min=" + String(d_min) + ", d_max=" + String(d_max) + " and d_iter=" + String(d_iter), "printlog.txt");
    end while;
    // Modelica.Utilities.Streams.print("setState_pT required " + String(iter) + " iterations to find d=" + String(d_iter) + " for input T=" + String(T) + " and p=" + String(p), "printlog.txt");
    assert(iter<iter_max, "setState_pTX did not converge, input was p=" + String(p) + " and T=" + String(T));

    state.p := p;
    state.T := T;
    state.d := d_iter;
    f := EoS.setHelmholtzDerivsFirst(d=state.d, T=state.T);
    state.s := EoS.s(f=f);
    state.h := EoS.h(f=f);
    state.u := EoS.u(f=f);
  end setState_pTX;


  redeclare function extends setState_phX
  "Return thermodynamic state as function of (p, h)"

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau "inverse reduced temperature";

    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    SpecificEnthalpy h_crit=fluidConstants[1].HCRIT0;

    EoS.HelmholtzDerivs f;
    SaturationProperties sat;
    MassFraction x "vapour quality";

    Density d_min;
    Density d_max;
    Density d_iter;
    Density d_iter_old;
    Temperature T_min;
    Temperature T_max;
    Temperature T_iter;
    Temperature T_iter_old;

    Real RES[2] "residual function vector";
    Real RSS "residual sum of squares";
    Real RSS_old "residual sum of squares";
    Real Jacobian[2,2] "Jacobian matrix";
    Real NS[2] "Newton step vector";
    Real grad[2] "gradient vector";
    Real slope;

    constant Real tolerance=1e-9 "tolerance for RSS";
    Integer iter = 0;
    constant Integer iter_max = 200;
    Real lambda(min=1e-3,max=1) = 1 "convergence speed, default=1";

    Boolean useBacktrack=false;
    Integer backtrack = 0;
    Real RSS_bt;
    Real lambda_bt;
    constant Real lambda_min = 0.1 "minimum for convergence speed";
    Real lambda_temp = 1 "temporary variable for convergence speed";
    constant Real alpha(min=0,max=1)=1e-4;
    Real rhs1;
    Real rhs2;
    Real a;
    Real b;
    Real Discriminant;

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(p <= p_crit, "setState_phX_error: pressure is higher than critical pressure");
      sat := setSat_p(p=p);
      assert(h >= sat.liq.h, "setState_phX_error: enthalpy is lower than saturated liquid enthalpy: this is single phase liquid");
      assert(h <= sat.vap.h, "setState_phX_error: enthalpy is higher than saturated vapor enthalpy: this is single phase vapor");
    else
      if (p < p_crit) then
        // two-phase possible, check region
        if (p>0.98*p_crit) or (p<300*p_trip) then
          // Modelica.Utilities.Streams.print("close to critical or triple point, get saturation properties from EoS", "printlog.txt");
          sat := setSat_p(p=p);
        else
          // do a simple check first, quite often this is sufficient
          sat.Tsat := Ancillary.saturationTemperature_p(p=p);

          sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
          f := EoS.setHelmholtzDerivsFirst(d=sat.liq.d, T=sat.Tsat, phase=1);
          sat.liq.h := EoS.h(f);

          sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
          f := EoS.setHelmholtzDerivsFirst(d=sat.vap.d, T=sat.Tsat, phase=1);
          sat.vap.h := EoS.h(f);

          if ((h > sat.liq.h - abs(0.02*sat.liq.h)) and (h < sat.vap.h + abs(0.02*sat.vap.h))) then
            // Modelica.Utilities.Streams.print("two-phase state or close to it, get saturation properties from EoS", "printlog.txt");
            sat := setSat_p(p=p);
          end if;
        end if;

        // Modelica.Utilities.Streams.print("phase boundary determined; sat.liq.h=" + String(sat.liq.h) + " and sat.vap.h=" + String(sat.vap.h), "printlog.txt");
        if (h < sat.liq.h) then
          // Modelica.Utilities.Streams.print("single phase liquid", "printlog.txt");
          state.phase := 1;

          d_min := sat.liq.d*0.98;
          d_iter:= sat.liq.d*1.02;
          d_max := fluidLimits.DMAX*1.10;

          T_min := fluidLimits.TMIN*0.99;
          T_max := sat.Tsat*1.02;
          T_iter:= sat.Tsat*0.98;
        elseif (h > sat.vap.h) then
          // Modelica.Utilities.Streams.print("single phase vapor", "printlog.txt");
          state.phase := 1;

          d_min := fluidLimits.DMIN;
          d_iter:= sat.vap.d/10;
          d_max := sat.vap.d*1.02;

          T_min := sat.Tsat*0.98;
          T_max := fluidLimits.TMAX*1.10;
          T_iter:= sat.Tsat*1.02;
        else
          // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
          state.phase := 2;
        end if;

      else
        state.phase := 1;
        // Modelica.Utilities.Streams.print("p>=p_crit, only single-phase possible", "printlog.txt");
        if (h<=h_crit) then
          // Modelica.Utilities.Streams.print("h<=h_crit, single-phase super-critical liquid-like region", "printlog.txt");
          d_min := d_crit*0.98;
          d_iter:= fluidLimits.DMAX*0.9;
          d_max := fluidLimits.DMAX*1.1;

          T_min := fluidLimits.TMIN;
          T_iter:= Ancillary.saturationTemperature_h_liq(h=h);
          T_max := T_iter*1.2;
        else
          // Modelica.Utilities.Streams.print("h>h_crit, single-phase super-critical vapour-like region", "printlog.txt");
          // due to the curvature, Newton will converge better when starting from the ideal gas region (low d, high T)
          d_min := fluidLimits.DMIN;
          d_iter:= d_crit/10;
          d_max := fluidLimits.DMAX*1.1;

          T_min := fluidLimits.TMIN;
          T_iter:= T_crit*1.3;
          T_max := fluidLimits.TMAX*1.1;
        end if;
      end if;
    end if;

    // phase and region determination finished !

    if (state.phase == 2) then
      // Modelica.Utilities.Streams.print("two-phase, SaturationProperties are already known", "printlog.txt");
      state.p := p;
      state.h := h;
      state.T := sat.Tsat;
      x := (h - sat.liq.h)/(sat.vap.h - sat.liq.h);
      state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
    else
      // Modelica.Utilities.Streams.print("single-phase, use 2D Newton-Raphson, start with d_iter=" + String(d_iter) + " (d_min=" + String(d_min) + " and d_max=" + String(d_max) + ") and T_iter=" + String(T_iter)+ " (T_min=" + String(T_min) + " and T_max=" + String(T_max) + ")", "printlog.txt");
      f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
      RES := {EoS.p(f)-p, EoS.h(f)-h};
      RSS := RES*RES/2;

      while ((RSS>tolerance) and (iter<iter_max)) loop
        iter := iter+1;

        // calculate Jacobian matrix, Newton Step vector, gradient vector and slope
        Jacobian := [EoS.dpdT(f), EoS.dpTd(f);
                     EoS.dhdT(f), EoS.dhTd(f)];
        NS := -Modelica.Math.Matrices.solve(Jacobian,RES);
        grad := RES*Jacobian;
        slope := grad*NS;
        assert(slope<0,"roudoff problem");

        // store old d_iter, T_iter and RSS
        d_iter_old := d_iter;
        T_iter_old := T_iter;
        RSS_old := RSS;

        // calculate new d_iter and T_iter using full Newton step (lambda=1)
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];

        // check bounds
        d_iter := max(d_iter, d_min);
        d_iter := min(d_iter, d_max);
        T_iter := max(T_iter, T_min);
        T_iter := min(T_iter, T_max);

        // calculate new residual vector and residual sum of squares
        f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
        RES := {EoS.p(f)-p, EoS.h(f)-h};
        RSS := RES*RES/2;
        // Modelica.Utilities.Streams.print("iter=" + String(iter) + " d_iter=" + String(d_iter) + " T_iter=" + String(T_iter) + " RES_p=" + String(RES[1]) + " RES_h=" + String(RES[2]) + " RSS=" + String(RSS), "printlog.txt");

        // if RSS is not decreasing fast enough, the full Newton step is not used
        // instead, the backtracking / linesearching loop tries to find lambda such that RSS decreases
        while useBacktrack and (backtrack<iter_max) and (lambda>=lambda_min) and not (RSS<=(RSS_old+alpha*lambda*slope)) loop
          backtrack := backtrack+1;

          // decrease lambda
          if (lambda==1) then
            lambda_temp := -slope/(2*(RSS-RSS_old-slope));
          else
            rhs1 := RSS   -RSS_old-lambda   *slope;
            rhs2 := RSS_bt-RSS_old-lambda_bt*slope;
            a := (           rhs1/lambda^2 -        rhs2/lambda_bt^2) / (lambda-lambda_bt);
            b := (-lambda_bt*rhs1/lambda^2 + lambda*rhs2/lambda_bt^2) / (lambda-lambda_bt);
            if (a <=0) then
              lambda_temp := -slope/(2*b);
            else
              Discriminant := b*b-3*a*slope;
              if (Discriminant<0) then
                lambda_temp := 0.5*lambda;
              elseif (b<=0) then
                lambda_temp := (-b+sqrt(Discriminant))/(3*a);
              else
                lambda_temp := -slope/(b+sqrt(Discriminant));
              end if;
            end if;
            // new lambda should be less or equal 0.5*previous lambda
            lambda_temp := if (lambda_temp>0.5*lambda) then 0.5*lambda else lambda_temp;
          end if;
          // store values for subsequent backtrack attempt
          lambda_bt := lambda;
          RSS_bt := RSS;

          // new lambda should be greater or equal 0.1*previous lambda
          lambda := max({lambda_temp, 0.1*lambda});

          d_iter := d_iter_old +lambda*NS[1];
          T_iter := T_iter_old +lambda*NS[2];
          // check bounds
          d_iter := max(d_iter, d_min);
          d_iter := min(d_iter, d_max);
          T_iter := max(T_iter, T_min);
          T_iter := min(T_iter, T_max);

          // calculate new residual vector and residual sum of squares
          f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
          RES := {EoS.p(f)-p, EoS.h(f)-h};
          RSS := RES*RES/2;
          // Modelica.Utilities.Streams.print("    backtrack attempt " + String(backtrack) + ": lambda= " + String(lambda) + " d_iter=" + String(d_iter) + " T_iter= " + String(T_iter) + " RSS= " + String(RSS), "printlog.txt");

        end while;
        // reset backtrack and lambda
        backtrack := 0;
        lambda_temp := 1;
        lambda := 1;

      end while;
      // Modelica.Utilities.Streams.print("setState_ph total iteration steps " + String(iter), "printlog.txt");
      assert(iter<iter_max, "setState_phX did not converge, input was p=" + String(p) + " and h=" + String(h));

      state.p := p;
      state.h := h;
      state.d := d_iter;
      state.T := T_iter;
      state.u := EoS.u(f);
      state.s := EoS.s(f);
    end if;

  end setState_phX;


  redeclare function extends setState_psX
  "Return thermodynamic state as function of (p, s)"

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau "inverse reduced temperature";

    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    SpecificEntropy s_crit=fluidConstants[1].SCRIT0;

    EoS.HelmholtzDerivs f;
    SaturationProperties sat;
    MassFraction x "vapour quality";

    Density d_min;
    Density d_max;
    Density d_iter;
    Density d_iter_old;
    Temperature T_min;
    Temperature T_max;
    Temperature T_iter;
    Temperature T_iter_old;

    Real RES[2] "residual function vector";
    Real RSS "residual sum of squares";
    Real RSS_old "residual sum of squares";
    Real Jacobian[2,2] "Jacobian matrix";
    Real NS[2] "Newton step vector";
    Real grad[2] "gradient vector";

    Real lambda(min=0.1,max=1) = 1 "convergence speed, default=1";
    constant Real tolerance=1e-9 "tolerance for RSS";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(p <= p_crit, "setState_psX_error: pressure is higher than critical pressure");
      sat := setSat_p(p=p);
      assert(s >= sat.liq.s, "setState_psX_error: entropy is lower than saturated liquid entropy: this is single phase liquid");
      assert(s <= sat.vap.s, "setState_psX_error: entropy is higher than saturated vapor entropy: this is single phase vapor");
    else
      if (p <= p_crit) then
        // two-phase possible, check region
        if (p>0.98*p_crit) or (p<300*p_trip) then
          // Modelica.Utilities.Streams.print("close to critical or triple point, get saturation properties from EoS", "printlog.txt");
          sat := setSat_p(p=p);
        else
          // do a simple check first, quite often this is sufficient
          sat.Tsat := Ancillary.saturationTemperature_p(p=p);

          sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
          f := EoS.setHelmholtzDerivsFirst(d=sat.liq.d, T=sat.Tsat, phase=1);
          sat.liq.s := EoS.s(f);

          sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
          f := EoS.setHelmholtzDerivsFirst(d=sat.vap.d, T=sat.Tsat, phase=1);
          sat.vap.s := EoS.s(f);

          if ((s > sat.liq.s - abs(0.05*sat.liq.s)) and (s < sat.vap.s + abs(0.05*sat.vap.s))) then
            // Modelica.Utilities.Streams.print("two-phase state or close to it, get saturation properties from EoS", "printlog.txt");
            sat := setSat_p(p=p);
          end if;
        end if;

        // Modelica.Utilities.Streams.print("phase boundary determined; sat.liq.s=" + String(sat.liq.s) + " and sat.vap.s=" + String(sat.vap.s), "printlog.txt");
        if (s < sat.liq.s) then
          // Modelica.Utilities.Streams.print("single phase liquid", "printlog.txt");
          state.phase := 1;
          d_min := sat.liq.d;
          d_max := 1.1*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
          d_iter:= sat.liq.d;
          T_min := fluidLimits.TMIN;
          T_max := sat.Tsat;
          T_iter:= sat.Tsat;
        elseif (s > sat.vap.s) then
          // Modelica.Utilities.Streams.print("single phase vapor", "printlog.txt");
          state.phase := 1;
          d_min := fluidLimits.DMIN;
          d_max := sat.vap.d;
          d_iter:= sat.vap.d/10;
          T_min := sat.Tsat;
          T_max := fluidLimits.TMAX;
          T_iter:= sat.Tsat;
        else
          // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
          state.phase := 2;
        end if;

      else
        state.phase := 1;
        // Modelica.Utilities.Streams.print("p>=p_crit, only single-phase possible", "printlog.txt");
        if (s<=s_crit) then
          // Modelica.Utilities.Streams.print("s<=s_crit, single-phase super-critical liquid-like region", "printlog.txt");
          d_min := d_crit;
          d_max := 1.1*fluidLimits.DMAX;
          d_iter:= 0.9*fluidLimits.DMAX;
          T_min := fluidLimits.TMIN;
          T_iter:= Ancillary.saturationTemperature_s_liq(s=s);
          T_max := 2*T_iter;
        else
          // Modelica.Utilities.Streams.print("s>s_crit, single-phase super-critical vapour-like region", "printlog.txt");
          d_min := fluidLimits.DMIN;
          d_max := fluidLimits.DMAX;
          d_iter:= p/(R*T_crit);
          T_min := T_crit;
          T_max := fluidLimits.TMAX;
          T_iter:= p/(R*d_crit);
        end if;
      end if;
    end if;

    // phase and region determination finished !

    if (state.phase == 2) then
      // force two-phase, SaturationProperties are already known
      state.p := p;
      state.s := s;
      state.T := sat.Tsat;
      x := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
      state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    else
      // single-phase, use 2D Newton-Raphson with linesearch and backtrack
      f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
      RES := {EoS.p(f)-p, EoS.s(f)-s};
      RSS := RES*RES/2;

      while ((RSS>tolerance) and (iter<iter_max)) loop
        iter := iter+1;

        // set up Jacobian matrix
        Jacobian := [EoS.dpdT(f), EoS.dpTd(f);
                     EoS.dsdT(f), EoS.dsTd(f)];

        // calculate vector of full Newton steps
        NS := -Modelica.Math.Matrices.solve(Jacobian,RES);

        // store old d_iter, T_iter and RSS
        d_iter_old := d_iter;
        T_iter_old := T_iter;
        RSS_old := RSS;

        // calculate new d_iter and T_iter using full Newton step
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];

        // check bounds
        d_iter := max(d_iter, 0.98*d_min);
        d_iter := min(d_iter, 1.02*d_max);
        T_iter := max(T_iter, 0.98*T_min);
        T_iter := min(T_iter, 1.02*T_max);

        // calculate new residual vector and residual sum of squares
        f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
        RES := {EoS.p(f)-p, EoS.s(f)-s};
        RSS := RES*RES/2;
      end while;
      // Modelica.Utilities.Streams.print("setState_ps total iteration steps " + String(iter), "printlog.txt");
      assert(iter<iter_max, "setState_psX did not converge, input was p=" + String(p) + " and s=" + String(s));

      state.p := p;
      state.s := s;
      state.d := d_iter;
      state.T := T_iter;
      state.u := EoS.u(f);
      state.h := EoS.h(f);
    end if;

  end setState_psX;


  redeclare function extends temperature
  "returns temperature from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output T
  algorithm
    T := state.T;
  end temperature;


  redeclare function extends density
  "returns density from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output d
  algorithm
    d := state.d;
  end density;


  redeclare function extends pressure
  "returns pressure from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output p
  algorithm
    p := state.p;
  end pressure;


  redeclare function extends specificEntropy
  "returns specificEntropy from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output h
  algorithm
    s := state.s;
  end specificEntropy;


  redeclare function extends specificEnthalpy
  "returns specificEnthalpy from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output h
  algorithm
    h := state.h;
  end specificEnthalpy;


  redeclare function vapourQuality "returns the vapour quality"
    // redeclare with algorithm based on d and T
    // previously only input state and output x were defined
    // optional input for saturation properties can save some time

    input ThermodynamicState state;
    input SaturationProperties sat=setSat_T(state.T);
    output MassFraction x;

protected
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Temperature T_crit=fluidConstants[1].criticalTemperature;

  algorithm
    // assert(state.T >= T_trip, "vapourQuality error: Temperature is lower than triple-point temperature");
    // assert(state.T <= T_crit, "vapourQuality error: Temperature is higher than critical temperature");

    if state.d <= sat.vap.d then
      x := 1;
    elseif state.d >= sat.liq.d then
      x := 0;
    else
      x := (1.0/state.d - 1.0/sat.liq.d)/(1.0/sat.vap.d - 1.0/sat.liq.d);
    end if;

  end vapourQuality;


  redeclare function extends specificHeatCapacityCp
  "returns the isobaric specific heat capcacity"
  //input state
  //output cp

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if (state.phase == 1) then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      cp := EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
    elseif (state.phase == 2) then
      // assert(false, "specificHeatCapacityCp warning: property not defined in two-phase region");
      cp := Modelica.Constants.inf; // division by zero
    end if;

  end specificHeatCapacityCp;


  redeclare function extends specificHeatCapacityCv
  "returns the isochoric specific heat capcacity"
  //input state
  //output cv

protected
    EoS.HelmholtzDerivs f;

    SaturationProperties sat;
    MassFraction x "vapour quality";
    DerPressureByTemperature dpT;
    EoS.HelmholtzDerivs fl;
    EoS.HelmholtzDerivs fv;

    DerEnergyByTemperature duT_liq;
    DerEnergyByTemperature duT_vap;
    DerDensityByTemperature ddT_liq;
    DerDensityByTemperature ddT_vap;
    DerVolumeByTemperature dvT_liq;
    DerVolumeByTemperature dvT_vap;
    DerFractionByTemperature dxTv;

  algorithm
    if (state.phase == 1) then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      cv := EoS.duTd(f);

    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T);
      x := (1/state.d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
      dpT := (sat.vap.s-sat.liq.s)/(1.0/sat.vap.d-1.0/sat.liq.d);

      fl := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.liq.d, phase=1);
      fv := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.vap.d, phase=1);

      duT_liq := EoS.duTd(fl)-EoS.dudT(fl)*EoS.dpTd(fl)/EoS.dpdT(fl) + EoS.dudT(fl)/EoS.dpdT(fl)*dpT;
      duT_vap := EoS.duTd(fv)-EoS.dudT(fv)*EoS.dpTd(fv)/EoS.dpdT(fv) + EoS.dudT(fv)/EoS.dpdT(fv)*dpT;
      ddT_liq := -EoS.dpTd(fl)/EoS.dpdT(fl) + 1.0/EoS.dpdT(fl)*dpT;
      ddT_vap := -EoS.dpTd(fv)/EoS.dpdT(fv) + 1.0/EoS.dpdT(fv)*dpT;
      dvT_liq := -1/sat.liq.d^2 * ddT_liq;
      dvT_vap := -1/sat.vap.d^2 * ddT_vap;
      dxTv :=(x*dvT_vap + (1 - x)*dvT_liq)/(1/sat.liq.d - 1/sat.vap.d);

      cv := duT_liq + dxTv*(sat.vap.u-sat.liq.u) + x*(duT_vap-duT_liq);
    end if;

  end specificHeatCapacityCv;


  redeclare function extends velocityOfSound
  "returns the speed or velocity of sound"
  //input state and
  //output a are inherited from PartialMedium

protected
    EoS.HelmholtzDerivs f;

  algorithm
    assert(state.phase <> 2, "velocityOfSound error: property not defined in two-phase region");
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
    a := sqrt(EoS.dpdT(f)-EoS.dpTd(f)*EoS.dsdT(f)/EoS.dsTd(f));
  end velocityOfSound;


  redeclare function extends isobaricExpansionCoefficient
  "returns 1/v*(dv/dT)@p=const"
  //input state
  //output beta

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if (state.phase == 1) then
      f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      beta := 1.0/state.d*EoS.dpTd(f)/EoS.dpdT(f);
    elseif (state.phase == 2) then
      beta := Modelica.Constants.small; // zero
    end if;
  end isobaricExpansionCoefficient;


  redeclare function extends isothermalCompressibility
  "returns -1/v*(dv/dp)@T=const"
  //input state and
  //output kappa are inherited from PartialMedium

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if (state.phase == 1) then
      f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      kappa := 1.0/state.d/EoS.dpdT(f);
    elseif (state.phase == 2) then
      kappa := Modelica.Constants.inf; // divide by zero
    end if;
  end isothermalCompressibility;


  redeclare function extends isentropicExponent "returns cp/cv"
  //input state and
  //output gamma are inherited from PartialMedium

protected
    EoS.HelmholtzDerivs f;
    SpecificHeatCapacity cp;
    SpecificHeatCapacity cv;

  algorithm
    if (state.phase == 1) then
      f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      cp := EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
      cv := EoS.duTd(f);
      gamma := cp/cv;
    elseif (state.phase == 2) then
      gamma := Modelica.Constants.inf;
    end if;
  end isentropicExponent;


  redeclare replaceable function extends dynamicViscosity
  "Returns dynamic Viscosity"
    // inherits input state and output eta

protected
    constant Real micro=1e-6;

  algorithm
    // assert(state.phase <> 2, "dynamicViscosity error: property not defined in two-phase region");

    // RefProp results are in µPa·s where µ means micro or 1E-6 but SI default is Pa·s
    eta := micro*(Transport.dynamicViscosity_dilute(state)
                + Transport.dynamicViscosity_initial(state)
                + Transport.dynamicViscosity_residual(state));

    annotation (Documentation(info="<html>
<p>
This model is identical to the RefProp VS1 or VS2 model.

The viscosity is split into three contributions: 
zero density (dilute gas) viscosity eta_0, 
initial density contribution eta_1
and residual contribution eta_r.

This allows to develop functions for each contribution seperately.
The so called background viscosity is the sum of initial and residual viscosity.

At the critical point and a small region around the critical point, the viscosity is enhanced. 
As this critical enhancement is small, it is neglected here.

Special thanks go to Eric W. Lemmon for answering all my emails 
and programming a special version of RefProp that outputs also intermediate values.

</p>

<dl>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
<dt>Vogel, E.; Küchenmeister, C. and Birch, E.</dt>
<dd> <b>Reference correlation of the viscosity of propane</b>.<br>
     Journal of Thermophysics (1998) 10, 417-426.<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
  end dynamicViscosity;


  redeclare replaceable function extends thermalConductivity
  "Return thermal conductivity"
    // inherits input state and output lambda

protected
    constant Real milli=1e-3;

  algorithm
    // assert(state.phase <> 2, "thermalConductivity error: property not defined in two-phase region");

    // RefProp results are in mW/m·K but SI default is W/m·K
    lambda := milli*(Transport.thermalConductivity_dilute(state)
                   + Transport.thermalConductivity_residual(state)
                   + Transport.thermalConductivity_critical(state));

    annotation (Documentation(info="<html>
  <p>
The thermal conductivity (TC) is split into three parts: ideal gas TC lamda_0, residual TC lambda_r and critical TC enhancement lambda_c.
Sometimes the residual TC is again split into two parts.
This allows to develop functions for each contribution seperately.
The sum of ideal gas TC and residual TC is called background TC.
Ideal gas TC depends on Temperature only and can be modelled by a quadratic function.
Residual TC is also modeled by a polynominal.
At the critical point TC becomes infinite; TC is enhanced for a large region around the critical point.
The critical enhancement can be described by various alternative approaches.
Here, the simplified approach as suggested by Olchowy and Sengers is implemented.

Special thanks go to Eric W. Lemmon for answering all my emails 
and programming a special version of RefProp that outputs also intermediate values.

</p>

</html>"));
  end thermalConductivity;


  redeclare replaceable function extends surfaceTension
  "Return surface tension sigma in the two phase region"
      // inherits input saturationProperties sat and output SurfaceTension sigma
      // this algorithm uses T only
      // liquid and vapour density are used in some mixture models

protected
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Temperature T_crit=fluidConstants[1].criticalTemperature;

    Real[size(surfaceTensionCoefficients.coeffs, 1)] a=
        surfaceTensionCoefficients.coeffs[:, 1];
    Real[size(surfaceTensionCoefficients.coeffs, 1)] n=
        surfaceTensionCoefficients.coeffs[:, 2];
    Real X "reduced temperature difference";

  algorithm
    assert(sat.Tsat >= T_trip, "vapourQuality error: Temperature is lower than triple-point temperature");
    assert(sat.Tsat <= T_crit, "vapourQuality error: Temperature is higher than critical temperature");

    X := (T_crit - sat.Tsat)/T_crit;
    sigma := sum(a[i]*X^n[i] for i in 1:size(a, 1));

    annotation (Documentation(info="<html>
  <p>
This is an implementation of the model as suggested by Somayajulu, G.R., 
which is an extension of the van der Waals surface tension correlation. 
The extended version has up to three terms with two parameters each.
</p>
<dl>
<dt>Somayajulu, G.R.</dt>
<dd> <b>A generalized equation for surface tension from the triple point to the critical point</b>.<br>
     International Journal of Thermophysics (1988) 9, 559-566.<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF00503154\">10.1007/BF00503154</a>
</dd>
<dt>Van der Waals, J.D.</dt>
<dd> <b>Thermodynamische Theorie der Kapillarität unter Voraussetzung stetiger Dichteänderung</b>.<br>
     Zeitschrift für Physikalische Chemie (1894) 13, 657-725.
</dd>
</dl>
</html>"));
  end surfaceTension;


  redeclare function extends bubbleEnthalpy
  "returns specificEnthalpy from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output hl
  algorithm
    hl := sat.liq.h;
  end bubbleEnthalpy;


  redeclare function extends dewEnthalpy
  "returns specificEnthalpy from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output hl
  algorithm
    hv := sat.vap.h;
  end dewEnthalpy;


  redeclare function extends dewDensity
  "returns density from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output hl
  algorithm
    dv := sat.vap.d;
  end dewDensity;


  redeclare function extends bubbleDensity
  "returns density from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output hl
  algorithm
    dl := sat.liq.d;
  end bubbleDensity;


  redeclare function extends pressure_dT
  // input, output and algorithm are inherited from PartialTwoPhaseMedium
  annotation (
    derivative=pressure_dT_der,
    inverse(d=density_pT(p=p, T=T, phase=phase),
            T=temperature_pd(p=p, d=d, phase=phase)));
  end pressure_dT;


  redeclare function extends density_pT
  // input, output and algorithm are inherited from PartialTwoPhaseMedium

  annotation (
    inverse(p=pressure_dT(d=d, T=T, phase=phase),
            T=temperature_pd(p=p, d=d, phase=phase)));
  end density_pT;


  redeclare function temperature_ps "returns temperature for given p and d"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Entropy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Temperature T "Temperature";

  algorithm
    T := temperature(setState_ps(p=p, s=s, phase=phase));

  annotation (
    inverse(p=pressure_Ts(T=T, s=s, phase=phase),
            s=specificEntropy_pT(p=p, T=T, phase=phase)));
  end temperature_ps;


  redeclare function specificEnthalpy_pT
  "returns specific enthalpy for given p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "specific enthalpy";

  algorithm
    h := specificEnthalpy(setState_pTX(p=p, T=T, phase=phase));

  annotation (
    inverse(T=temperature_ph(p=p, h=h, phase=phase)));
  end specificEnthalpy_pT;


  redeclare function temperature_ph "returns temperature for given p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Enthalpy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Temperature T "Temperature";

  algorithm
    T := temperature(setState_ph(p=p, h=h, phase=phase));

  annotation (
    inverse(h=specificEnthalpy_pT(p=p, T=T, phase=phase)));
  end temperature_ph;


  redeclare function specificEnthalpy_ps
  "returns specific enthalpy for a given p and s"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Entropy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "specific enthalpy";

  algorithm
    h := specificEnthalpy(setState_psX(p=p, s=s, phase=phase));

  annotation (
    inverse(s=specificEntropy_ph(p=p, h=h, phase=phase)));
  end specificEnthalpy_ps;


  redeclare function density_ph "returns density for given p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Enthalpy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Density d "Temperature";

  algorithm
    d := density(setState_ph(p=p, h=h, phase=phase));

  annotation (
    derivative=density_ph_der,
    inverse(h=specificEnthalpy_pd(p=p, d=d, phase=phase)));
  end density_ph;


  redeclare function extends density_derp_h
  "returns density derivative (dd/dp)@h=const"
  //input state
  //output ddph

protected
    EoS.HelmholtzDerivs f;

    SaturationProperties sat;
    MassFraction x "vapour quality";
    DerTemperatureByPressure dTp;
    EoS.HelmholtzDerivs fl;
    EoS.HelmholtzDerivs fv;

    DerEnthalpyByPressure dhp_liq;
    DerEnthalpyByPressure dhp_vap;
    DerDensityByPressure ddp_liq;
    DerDensityByPressure ddp_vap;
    DerVolumeByPressure dvp_liq;
    DerVolumeByPressure dvp_vap;
    DerVolumeByPressure dvph;
    DerFractionByPressure dxph;

  algorithm
    if (state.phase == 1) then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
      ddph := 1.0/(EoS.dpdT(f) - EoS.dpTd(f)*EoS.dhdT(f)/EoS.dhTd(f));
    elseif (state.phase == 2) then
      sat := setSat_T(T=state.T);
      x := (state.h - sat.liq.h)/(sat.vap.h - sat.liq.h);
      dTp := (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

      fl := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.liq.d, phase=1);
      fv := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.vap.d, phase=1);

      dhp_liq := EoS.dhdT(fl)/EoS.dpdT(fl) + (EoS.dhTd(fl)-EoS.dhdT(fl)*EoS.dpTd(fl)/EoS.dpdT(fl))*dTp;
      dhp_vap := EoS.dhdT(fv)/EoS.dpdT(fv) + (EoS.dhTd(fv)-EoS.dhdT(fv)*EoS.dpTd(fv)/EoS.dpdT(fv))*dTp;
      ddp_liq := 1.0/EoS.dpdT(fl) - EoS.dpTd(fl)/EoS.dpdT(fl)*dTp;
      ddp_vap := 1.0/EoS.dpdT(fv) - EoS.dpTd(fv)/EoS.dpdT(fv)*dTp;
      dvp_liq := -1.0/sat.liq.d^2 * ddp_liq;
      dvp_vap := -1.0/sat.vap.d^2 * ddp_vap;
      dxph :=(x*dhp_vap + (1 - x)*dhp_liq)/(sat.liq.h - sat.vap.h);

      dvph := dvp_liq + dxph*(1.0/sat.vap.d-1.0/sat.liq.d) + x*(dvp_vap-dvp_liq);
      ddph := -state.d*state.d*dvph;
    end if;
  end density_derp_h;


  redeclare function extends density_derh_p
  "returns density derivative (dd/dh)@p=const"
  //input state
  //output ddhp

protected
    EoS.HelmholtzDerivs f;
    SaturationProperties sat;

  algorithm
    if (state.phase == 1) then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
      ddhp := 1.0/(EoS.dhdT(f) - EoS.dhTd(f)*EoS.dpdT(f)/EoS.dpTd(f));
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T);
      // dvhp = (v"-v')/(h"-h')
      // ddhp = -d^2 * dvhp
      ddhp := -state.d^2*(1/sat.liq.d-1/sat.vap.d)/(sat.liq.h-sat.vap.h);
    end if;
  end density_derh_p;


  redeclare function extends density_derp_T
  "returns density derivative (dd/dp)@T=const"
  //input state and output ddpT are inherited

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if (state.phase == 1) then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      ddpT := 1.0/EoS.dpdT(f);
    elseif (state.phase == 2) then
      ddpT := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derp_T;


  redeclare function extends density_derT_p
  "returns density derivative (dd/dT)@p=const"
  //input state and output ddTp are inherited

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if (state.phase == 1) then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      ddTp := -EoS.dpTd(f)/EoS.dpdT(f);
    elseif (state.phase == 2) then
      ddTp := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derT_p;


  redeclare function extends saturationTemperature

protected
    SaturationProperties sat = setSat_p(p=p);

  algorithm
    T := sat.Tsat;
  end saturationTemperature;


  redeclare function saturationTemperature_derp "returns (dT/dp)@sat"
  // does not extend, because base class output has wrong units
  input AbsolutePressure p;
  output DerTemperatureByPressure dTp;

protected
   SaturationProperties sat=setSat_p(p=p);

  algorithm
    // Clausius-Clapeyron equation
    dTp := (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);
  end saturationTemperature_derp;


  redeclare function extends saturationPressure

protected
    SaturationProperties sat = setSat_T(T=T);

  algorithm
    p := sat.psat;
  end saturationPressure;


  redeclare function extends dBubbleDensity_dPressure
  "Return bubble point density derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output ddldp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.liq.T);
    DerDensityByPressure ddpT = 1.0/EoS.dpdT(f);
    DerDensityByTemperature ddTp = -EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    ddldp := ddpT + ddTp*dTp;
  end dBubbleDensity_dPressure;


  redeclare function extends dDewDensity_dPressure
  "Return dew point density derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output ddvdp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.vap.T);
    DerDensityByPressure ddpT = 1.0/EoS.dpdT(f);
    DerDensityByTemperature ddTp = -EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    ddvdp := ddpT + ddTp*dTp;
  end dDewDensity_dPressure;


  redeclare function extends dBubbleEnthalpy_dPressure
  "Return bubble point enthalpy derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dhldp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.liq.T);
    DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
    DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    dhldp := dhpT + dhTp*dTp;
  end dBubbleEnthalpy_dPressure;


  redeclare function extends dDewEnthalpy_dPressure
  "Return dew point enthalpy derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dhvdp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.vap.T);
    DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
    DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    dhvdp := dhpT + dhTp*dTp;
  end dDewEnthalpy_dPressure;


  redeclare function extends specificEnthalpy_dT
  // input, output and algorithm are inherited from PartialTwoPhaseMedium
  annotation(derivative=specificEnthalpy_dT_der);
  end specificEnthalpy_dT;
end PartialHelmholtzMedium;
