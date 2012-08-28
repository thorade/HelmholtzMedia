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
    Density d "Density of medium";
    Temperature T "Temperature of medium";
    AbsolutePressure p "Absolute pressure of medium";
    SpecificEnthalpy h "Specific enthalpy of medium";
    SpecificEnergy u "Specific inner energy of medium";
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

    // Boolean verbose=true;
protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    EoS.HelmholtzDerivs fl;
    EoS.HelmholtzDerivs fv;

    Real RES_pl;
    Real f1dx;
    //Real f1dy;
    Real f1dz;

    Real RES_pv;
    //Real f2dx;
    Real f2dy;
    Real f2dz;

    Real RES_g;
    Real f3dx;
    Real f3dy;
    Real f3dz;

    Real det "determinant of Jacobi matrix";
    Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
    Real tolerance=1e-3
    "tolerance for sum of RES_pl (in Pa), RES_pv (in Pa) and RES_g (in J/kg)";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
  if ((p>=p_trip) and (p<p_crit)) then
    // calculate start values
    // sat.Tsat  := 1/(1/T_crit - (1/T_trip-1/T_crit)/log(p_crit/p_trip)*log(p/p_crit));
    sat.Tsat := Ancillary.saturationTemperature_p(p=p);
    sat.liq.d := 1.02*Ancillary.bubbleDensity_T(T=sat.Tsat);
    sat.vap.d := 0.98*Ancillary.dewDensity_T(T=sat.Tsat);

    // calculate residuals
    fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, phase=1);
    fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, phase=1);
    RES_pl := fl.d*fl.T*fl.R*(1+fl.delta*fl.rd) - p;   // f1
    RES_pv := fv.d*fv.T*fv.R*(1+fv.delta*fv.rd) - p;   // f2
    RES_g  := fl.T*fl.R*((fl.i+fl.r)+(1+fl.delta*fl.rd))
            - fv.T*fv.R*((fv.i+fv.r)+(1+fv.delta*fv.rd));  // f3

    while (((abs(RES_pl) + abs(RES_pv) + abs(RES_g)) > tolerance) and (iter<iter_max)) loop
      iter := iter+1;

      // calculate gradients of f1,f2,f3 with respect to x,y,z
      f1dx := fl.T*fl.R*(1+2*fl.delta*fl.rd+fl.delta^2*fl.rdd);
      //f1dy := 0;
      f1dz := fl.d*fl.R*(1+fl.delta*fl.rd-fl.delta*fl.tau*fl.rtd);
      //f2dx := 0;
      f2dy := fv.T*fv.R*(1+2*fv.delta*fv.rd+fv.delta^2*fv.rdd);
      f2dz := fv.d*fl.R*(1+fv.delta*fv.rd-fv.delta*fv.tau*fv.rtd);
      f3dx := +fl.R*fl.T/fl.d*(1 + 2*fl.delta*fl.rd + fl.delta^2*fl.rdd);
      f3dy := -fv.R*fv.T/fv.d*(1 + 2*fv.delta*fv.rd + fv.delta^2*fv.rdd);
      f3dz := fl.R*(-fl.tau*(fl.it+fl.rt) + (fl.i+fl.r) + (1+fl.delta*fl.rd-fl.delta*fl.tau*fl.rtd))
            - fv.R*(-fv.tau*(fv.it+fv.rt) + (fv.i+fv.r) + (1+fv.delta*fv.rd-fv.delta*fv.tau*fv.rtd));

      // calculate determinant of Jacobi matrix
      det := f1dx*f2dy*f3dz -f3dx*f2dy*f1dz - f3dy*f2dz*f1dx;

      /* // print for debugging
    Modelica.Utilities.Streams.print(" ", "printlog.txt");
    Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
    Modelica.Utilities.Streams.print("sat.liq.d=" + String(sat.liq.d) + " and sat.vap.d=" + String(sat.vap.d) + " and sat.Tsat=" + String(sat.Tsat), "printlog.txt");
    Modelica.Utilities.Streams.print("RES_pl=" + String(RES_pl) + " and RES_pv=" + String(RES_pv) + " and RES_g=" + String(RES_g), "printlog.txt");
    Modelica.Utilities.Streams.print("f1dx=" + String(f1dx) + " and f1dz=" + String(f1dz), "printlog.txt");
    Modelica.Utilities.Streams.print("f2dy=" + String(f2dy) + " and f2dz=" + String(f2dz), "printlog.txt");
    Modelica.Utilities.Streams.print("f3dx=" + String(f3dx) + " and f3dy=" + String(f3dy), "printlog.txt");
    Modelica.Utilities.Streams.print("det(J)=" + String(det), "printlog.txt"); */

      // calculate better sat.liq.d, sat.vap.d and sat.Tsat
      sat.liq.d := sat.liq.d - gamma/det*(RES_pl*(f2dy*f3dz-f2dz*f3dy) +RES_pv*(f1dz*f3dy-0)         +RES_g*(0-f1dz*f2dy));
      sat.vap.d := sat.vap.d - gamma/det*(RES_pl*(f2dz*f3dx-0)         +RES_pv*(f1dx*f3dz-f1dz*f3dx) +RES_g*(0-f1dx*f2dz));
      sat.Tsat  := sat.Tsat  - gamma/det*(RES_pl*(0-f2dy*f3dx)         +RES_pv*(0-f1dx*f3dy)         +RES_g*(f1dx*f2dy-0));

      // check bounds
      sat.liq.d := max(sat.liq.d, d_crit);
      sat.liq.d := min(sat.liq.d, fluidLimits.DMAX);
      sat.vap.d := max(sat.vap.d, fluidLimits.DMIN);
      sat.vap.d := min(sat.vap.d, d_crit);
      sat.Tsat := max(sat.Tsat, T_trip);
      sat.Tsat := min(sat.Tsat, T_crit);

      // calculate new residual
      fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, phase=1);
      fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, phase=1);
      RES_pl := fl.d*fl.T*fl.R*(1+fl.delta*fl.rd) - p;   // f1
      RES_pv := fv.d*fv.T*fv.R*(1+fv.delta*fv.rd) - p;   // f2
      RES_g  := fl.T*fl.R*(fl.i+fl.r + fl.delta*fl.rd) - fv.T*fv.R*(fv.i+fv.r + fv.delta*fv.rd);  // f3
    end while;
    // if verbose then Modelica.Utilities.Streams.print("setSat_p total iteration steps " + String(iter), "printlog.txt"); end if;
    // Modelica.Utilities.Streams.print("setSat_p total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setSat_p did not converge, input was p=" + String(p));

    sat.psat := p;
    sat.liq := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);

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
  else
    // assert(p >= p_trip, "setSat_p error: pressure is lower than triple-point pressure");
    // p<p_trip: this does not make sense: if p is below the triple pressure, the medium is single phase vapour
    // anyway, during initialization (at time=0) p=0 may happen
    // fluidLimit values are returned
    sat.psat  := p;
    sat.Tsat  := T_trip;
    sat.liq.d := fluidLimits.DMAX;
    sat.vap.d := fluidLimits.DMIN;

  end if;
  end setSat_p;


  function setSat_d
  "iterative calculation of saturation properties from EoS with Newton-Raphson algorithm"
    input Density d;
    output SaturationProperties sat;

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real tau(unit="1") "inverse reduced temperature";

    Density d_min =  fluidLimits.DMIN;
    Density d_max =  fluidLimits.DMAX;

    EoS.HelmholtzDerivs fl;
    EoS.HelmholtzDerivs fv;

    Real RES_p;
    Real RES_g;
    Real dpdd;
    Real dpdT;
    Real dgdd;
    Real dgdT;
    Real det "determinant of Jacobi matrix";
    Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
    Real tolerance=1e-3
    "tolerance for sum of RES_p (in Pa) and RES_g (in J/kg)";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("setSat_d: d="+String(d),"printlog.txt");

    sat.Tsat := Ancillary.saturationTemperature_d(d=d);
    if (d<d_crit+tolerance) then
      // Modelica.Utilities.Streams.print("d<d_crit: input is on vapour side: find d_liq and T_sat", "printlog.txt");
      sat.liq.d := Ancillary.bubbleDensity_T(sat.Tsat);
      sat.vap.d := d; // d'' is a constant

      // calculate residuals: liq-vap
      fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, phase=1);
      fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, phase=1);
      RES_p := fl.d*fl.T*fl.R*(1+fl.delta*fl.rd)
             - fv.d*fv.T*fv.R*(1+fv.delta*fv.rd);
      RES_g := fl.T*fl.R*((fl.i+fl.r)+(1+fl.delta*fl.rd))
             - fv.T*fv.R*((fv.i+fv.r)+(1+fv.delta*fv.rd));

      while (abs(RES_p) + abs(RES_g) > tolerance) and iter<iter_max loop
        // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
        // Modelica.Utilities.Streams.print("sat.liq.d=" + String(sat.liq.d) + "  and dpdd=" + String(dpdd) + " and dgdd=" + String(dgdd), "printlog.txt");
        // Modelica.Utilities.Streams.print(" sat.Tsat=" + String(sat.Tsat)  + " and dpdT=" + String(dpdT) + " and dgdT=" + String(dgdT), "printlog.txt");
        iter := iter+1;

        // calculate gradients ragrding d_liq and T
        dpdd := fl.T*fl.R*(1+2*fl.delta*fl.rd + fl.delta^2*fl.rdd);
        dpdT := fl.d*fl.R*(1+fl.delta*fl.rd-fl.delta*fl.tau*fl.rtd)
              - fv.d*fv.R*(1+fv.delta*fv.rd-fv.delta*fv.tau*fv.rtd);
        dgdd := fl.T*fl.R/fl.d*(1+2*fl.delta*fl.rd + fl.delta^2*fl.rdd);
        dgdT := fl.R*(-fl.tau*(fl.it+fl.rt) +(fl.i+fl.r) +(1+fl.delta*fl.rd-fl.delta*fl.tau*fl.rtd))
              - fv.R*(-fv.tau*(fv.it+fv.rt) +(fv.i+fv.r) +(1+fv.delta*fv.rd-fv.delta*fv.tau*fv.rtd));

        // calculate determinant of Jacobi matrix det=ad-bc
        det := dpdd*dgdT-dpdT*dgdd;

        // calculate better values for sat.liq.d and sat.Tsat
        sat.liq.d := sat.liq.d -gamma/det*(+dgdT*RES_p -dpdT*RES_g);
        sat.Tsat  := sat.Tsat  -gamma/det*(-dgdd*RES_p +dpdd*RES_g);

        // check bounds
        sat.liq.d := max(sat.liq.d, d_crit);
        sat.liq.d := min(sat.liq.d, d_max);
        sat.Tsat  := max(sat.Tsat,  T_trip);
        sat.Tsat  := min(sat.Tsat,  T_crit);

        // calculate new residuals: liq-vap
        fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, phase=1);
        fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, phase=1);
        RES_p := fl.d*fl.T*fl.R*(1+fl.delta*fl.rd)
               - fv.d*fv.T*fv.R*(1+fv.delta*fv.rd);
        RES_g := fl.T*fl.R*((fl.i+fl.r)+(1+fl.delta*fl.rd))
               - fv.T*fv.R*((fv.i+fv.r)+(1+fv.delta*fv.rd));
      end while;

    elseif (d>d_crit-tolerance) then
      // Modelica.Utilities.Streams.print("d>d_crit: input is on liquid side: find d_vap and T_sat", "printlog.txt");
      sat.vap.d := Ancillary.dewDensity_T(sat.Tsat);
      sat.liq.d := d; // d' is a constant

      // calculate residuals: vap-liq
      fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, phase=1);
      fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, phase=1);
      RES_p := fv.d*fv.T*fv.R*(1+fv.delta*fv.rd)
             - fl.d*fl.T*fl.R*(1+fl.delta*fl.rd);
      RES_g := fv.T*fv.R*((fv.i+fv.r)+(1+fv.delta*fv.rd))
             - fl.T*fl.R*((fl.i+fl.r)+(1+fl.delta*fl.rd));

      while (abs(RES_p) + abs(RES_g) > tolerance) and iter<iter_max loop
        // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
        // Modelica.Utilities.Streams.print("sat.vap.d=" + String(sat.vap.d) + "  and dpdd=" + String(dpdd) + " and dgdd=" + String(dgdd), "printlog.txt");
        // Modelica.Utilities.Streams.print(" sat.Tsat=" + String(sat.Tsat)  + " and dpdT=" + String(dpdT) + " and dgdT=" + String(dgdT), "printlog.txt");
        iter := iter+1;

        // calculate gradients ragrding d_liq and T
        dpdd := fv.T*fv.R*(1+2*fv.delta*fv.rd + fv.delta^2*fv.rdd);
        dpdT := fv.d*fv.R*(1+fv.delta*fv.rd-fv.delta*fv.tau*fv.rtd)
              - fl.d*fl.R*(1+fl.delta*fl.rd-fl.delta*fl.tau*fl.rtd);
        dgdd := fv.T*fv.R/fv.d*(1+2*fv.delta*fv.rd + fv.delta^2*fv.rdd);
        dgdT := fv.R*(-fv.tau*(fv.it+fv.rt) +(fv.i+fv.r) +(1+fv.delta*fv.rd-fv.delta*fv.tau*fv.rtd))
              - fl.R*(-fl.tau*(fl.it+fl.rt) +(fl.i+fl.r) +(1+fl.delta*fl.rd-fl.delta*fl.tau*fl.rtd));

        // calculate determinant of Jacobi matrix det=ad-bc
        det := dpdd*dgdT-dpdT*dgdd;

        // calculate better values for sat.vap.d and sat.Tsat
        sat.vap.d := sat.vap.d -gamma/det*(+dgdT*RES_p -dpdT*RES_g);
        sat.Tsat  := sat.Tsat  -gamma/det*(-dgdd*RES_p +dpdd*RES_g);

        // check bounds
        sat.vap.d := max(sat.vap.d, d_min);
        sat.vap.d := min(sat.vap.d, d_crit);
        sat.Tsat  := max(sat.Tsat,  T_trip);
        sat.Tsat  := min(sat.Tsat,  T_crit);

        // calculate new residuals: liq-vap
        fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, phase=1);
        fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, phase=1);
        RES_p := fv.d*fv.T*fv.R*(1+fv.delta*fv.rd)
               - fl.d*fl.T*fl.R*(1+fl.delta*fl.rd);
        RES_g := fv.T*fv.R*((fv.i+fv.r)+(1+fv.delta*fv.rd))
               - fl.T*fl.R*((fl.i+fl.r)+(1+fl.delta*fl.rd));
      end while;

    else
      // Modelica.Utilities.Streams.print("d=d_crit: return critical Temperature", "printlog.txt");
      sat.Tsat  := T_crit;
      sat.liq.d := d_crit;
      sat.vap.d := d_crit;
    end if;
    // Modelica.Utilities.Streams.print("setSat_d total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setSat_d did not converge, input was d=" + String(d));

    sat.liq  := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap  := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);
    sat.psat := sat.liq.p;

  end setSat_d;


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
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Real delta(unit="1")=d/d_crit "reduced density";
    Real tau(unit="1")=T_crit/T "inverse reduced temperature";
    EoS.HelmholtzDerivs
                    f;

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
      x := (1/d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
      state.p := sat.psat;
      state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
    else
      // force single-phase
      f.i   := EoS.f_i(delta=delta, tau=tau);
      f.it  := EoS.f_it(delta=delta, tau=tau);
      f.r   := EoS.f_r(delta=delta, tau=tau);
      f.rt  := EoS.f_rt(delta=delta, tau=tau);
      f.rd  := EoS.f_rd(delta=delta, tau=tau);
      state.p := d*T*R*(1+delta*f.rd);
      state.h :=   T*R*(tau*(f.it + f.rt) + (1+delta*f.rd));
      state.u :=   T*R*(tau*(f.it + f.rt));
      state.s :=     R*(tau*(f.it + f.rt) - f.i - f.r);
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
    state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
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
    state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
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
    Real delta(unit="1") "reduced density";
    constant Real tau(unit="1")=T_crit/T "inverse reduced temperature";
    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    EoS.HelmholtzDerivs f;
    SaturationProperties sat;

    /* // RKS: Redlich-Kwong-Soave (see Soave 1979)
  constant Real omega = fluidConstants[1].acentricFactor;
  constant Real m = 0.480 + 1.574*omega - 0.176*omega^2;
  Real a = 0.42747*R^2*T_crit^2/p_crit*(1 + m*(1 - sqrt(T/T_crit)))^2;
  Real b = 0.08664*R*T_crit/p_crit;
  Real A = a*p/(R^2*T^2);
  Real B = b*p/(R*T);
  Real r = (A-B-B^2)-1/3;
  Real q = -2/27 + 1/3*(A-B-B^2) - A*B;
  Real D = (r/3)^3 + (q/2)^2 "discriminant";
  Real u;
  Real Y1;
  Real Y2;
  Real Y3;
  Real Theta;
  Real phi;
  import Modelica.Constants.pi;
  */

    Density d_min;
    Density d_max;
    Density d_med;
    Density d_iter;
    Real RES_p;
    Real RES_min;
    Real RES_max;
    Real RES_med;
    Real dpdd;
    constant Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
    constant Real tolerance=1e-9 "relativ tolerance for RES_p";
    Integer iter = 0;
    constant Integer iter_max=200;
    Boolean RiddersIsInitialized=false;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("p=" + String(p) + " and T=" + String(T), "printlog.txt");
    assert(phase <> 2, "setState_pTX_error: pressure and temperature are not independent variables in two-phase state");
    state.phase := 1;

    if (T < T_crit) then
      // determine p_sat
      sat.psat := Ancillary.saturationPressure_T(T=T);
      if (p > 1.02*sat.psat) then
        sat.liq.d := 0.97*Ancillary.bubbleDensity_T(T=T);
      elseif (p < 0.98*sat.psat) then
        sat.vap.d := 1.03*Ancillary.dewDensity_T(T=T);
      else
        // Modelica.Utilities.Streams.print("close to saturation boundary, get saturation properties from EoS", "printlog.txt");
        sat := setSat_T(T=T);
      end if;
      // Modelica.Utilities.Streams.print("sat.psat=" + String(sat.psat), "printlog.txt");

      // determine region
      if (p > sat.psat) then
        // Modelica.Utilities.Streams.print("single phase liquid: d is between dliq and rho_max", "printlog.txt");
        d_min  := sat.liq.d;
        d_max  := 1.1*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
        d_iter := sat.liq.d;
      elseif (p < sat.psat) then
        // Modelica.Utilities.Streams.print("single phase vapor: d is between 0 and dvap", "printlog.txt");
        d_min  := fluidLimits.DMIN;
        d_max  := sat.vap.d;
        d_iter := sat.vap.d/10;
      else
        // this should not happen
        assert(p <> sat.psat, "setState_pTX_error: pressure equals saturation pressure");
      end if;
    else
      // Modelica.Utilities.Streams.print("T>T_crit: d is between dmin and dmax", "printlog.txt");
      d_min  := fluidLimits.DMIN;
      d_max  := 1.1*fluidLimits.DMAX;
      d_iter := d_crit/50;
    end if;

    /* // get density start value from Redlich-Kwong-Soave (see Span 2000, section "3.3.1 Calculations based on pT" )  
  Modelica.Utilities.Streams.print("RKS discriminant D=" + String(D), "printlog.txt");
  if (D >= 0) then
    u := (sqrt(D)-(q/2))^(1/3);
    Y1 := u-r/(3*u);
    Modelica.Utilities.Streams.print("RKS has one root (Y1=" + String(Y1) + ")", "printlog.txt");
    d_iter := p/(R*T*(Y1+1/3));
  elseif ((abs(D) < 1e-8) and (abs(r) < 1e-3)) then
    // Modelica.Utilities.Streams.print("close to critical region, use critical density");
    d_iter := d_crit;
  else
    // case D<0
    Theta := sqrt(-r^3/27);
    phi := acos(-q/(2*Theta));
    Y1 := 2*Theta^(1/3)*cos(phi/3);
    Y2 := 2*Theta^(1/3)*cos(phi/3+2*pi/3);
    Y3 := 2*Theta^(1/3)*cos(phi/3+4*pi/3);
     Modelica.Utilities.Streams.print("RKS has three possible roots(Y1=" + String(Y1) + ", Y2=" + String(Y2) + ", Y3=" + String(Y3), "printlog.txt");
    if (T <= T_crit) then
      Modelica.Utilities.Streams.print("T<T_crit: multiple roots due to phase boundary", "printlog.txt");
      Modelica.Utilities.Streams.print("d(Y1)=" + String(p/(R*T*(Y1+1/3))) + ", d(Y2)=" + String(p/(R*T*(Y2+1/3))) + ", d(p/(R*T*(Y3+1/3)))=" + String(Y3), "printlog.txt");
      if (p > sat.psat) then
        Y1 := min(Y1,Y2);
        Y1 := min(Y1,Y3);
        d_iter := p/(R*T*(Y1+1/3));
      elseif (p < sat.psat) then
        Y1 := max(Y1,Y2);
        Y1 := max(Y1,Y3);
        d_iter := p/(R*T*(Y1+1/3));
      else
        // this should not happen
        assert(p <> sat.psat, "setState_pTX error: pressure equals saturation pressure");
      end if;
    else
      Modelica.Utilities.Streams.print("T>T_crit: multiple roots can occur, but two of the roots result in negative densities", "printlog.txt");
      Modelica.Utilities.Streams.print("d(Y1)=" + String(p/(R*T*(Y1+1/3))) + ", d(Y2)=" + String(p/(R*T*(Y2+1/3))) + ", d(p/(R*T*(Y3+1/3)))=" + String(Y3), "printlog.txt");
      d_iter := max(p/(R*T*(Y1+1/3)), p/(R*T*(Y2+1/3)));
      d_iter := max(p/(R*T*(Y3+1/3)), d_iter);
    end if;
  end if;
   Modelica.Utilities.Streams.print("RKS finished, d_iter=" + String(d_iter), "printlog.txt");
  // check bounds, RKS is not very accurate for VLE densities
  d_iter := max(d_min,d_iter);
  d_iter := min(d_max,d_iter);
  */

    // Modelica.Utilities.Streams.print("start Newton with d_min=" + String(d_min) + ", d_max=" + String(d_max) + " and d_iter=" + String(d_iter), "printlog.txt");
    // calculate RES_p
    delta := d_iter/d_crit;
    f.rd  := EoS.f_rd(delta=delta, tau=tau);
    RES_p := d_iter*T*R*(1+delta*f.rd) - p;

    while ((abs(RES_p/p) > tolerance) and (iter<iter_max)) loop
      iter := iter+1;

      // calculate gradient with respect to density
      f.rdd := EoS.f_rdd(delta=delta, tau=tau);
      dpdd := T*R*(1+2*delta*f.rd+delta^2*f.rdd);

      // print for Newton debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter) + ", current d_iter=" + String(d_iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("RES_p=" + String(RES_p) + " and dpdd=" + String(dpdd), "printlog.txt");

      // calculate better d_iter from Newton
      d_iter := d_iter - gamma/dpdd*RES_p;

      // check bounds, if out of bounds use Ridders
        if (d_iter<d_min) or (d_iter>d_max) then
          // Modelica.Utilities.Streams.print("d_iter out of bounds, fallback to Ridders' method, step=" + String(iter) + ", d_iter=" + String(d_iter), "printlog.txt");
          if not RiddersIsInitialized then
            // calculate RES_p for d_min
            delta := d_min/d_crit;
            f.rd  := EoS.f_rd(delta=delta, tau=tau);
            RES_min := d_min*T*R*(1+delta*f.rd) - p;
            // calculate RES_p for d_max
            delta := d_max/d_crit;
            f.rd  := EoS.f_rd(delta=delta, tau=tau);
            RES_max := d_max*T*R*(1+delta*f.rd) - p;
            // Modelica.Utilities.Streams.print("initialize Ridders, RES_min=" + String(RES_min) + " and RES_max=" + String( RES_max), "printlog.txt");
            RiddersIsInitialized := true;
          end if;
          if (RES_min*RES_max<0) then
            // calculate RES_p for d_med
            d_med := (d_max+1*d_min)/2;
            delta := d_med/d_crit;
            f.rd  := EoS.f_rd(delta=delta, tau=tau);
            RES_med := d_med*T*R*(1+delta*f.rd) - p;
            // find better d_iter by Ridders' method
            d_iter := d_med + (d_med-d_min)*sign(RES_min-RES_max)*RES_med/sqrt(RES_med^2-RES_min*RES_max);
            // calculate new RES_s
            delta := d_iter/d_crit;
            f.rd  := EoS.f_rd(delta=delta, tau=tau);
            RES_p := d_iter*T*R*(1+delta*f.rd) - p;
            // thighten the bounds
            if (RES_p*RES_med<=0) then
              // opposite sign, d_med and d_iter bracket the root
              d_min := d_iter;
              RES_min := RES_p;
              d_max := d_med;
              RES_max := RES_med;
            else
              if (RES_p*RES_min<0) then
                d_max := d_iter;
                RES_max := RES_p;
              elseif (RES_p*RES_max<0) then
                d_min := d_iter;
                RES_min := RES_p;
              else
                assert(false,"setState_pTX: this should never happen");
              end if;
            end if;
          // Modelica.Utilities.Streams.print("Ridders' method: new brackets d_min=" + String(d_min) + ", d_max=" + String(d_max), "printlog.txt");
          else
            if (abs(RES_min/p)<tolerance) then
              d_iter:= d_min;
            elseif (abs(RES_max/p)<tolerance) then
              d_iter:=d_max;
            else
              assert(false, "setState_pTX: d_min=" +String(d_min) + " and d_max =" + String(d_max) + " did not bracket the root, input was p=" + String(p) + " and T=" + String(T));
            end if;
          end if;
        else
          // d_iter from Newton is within bounds
          // calculate new RES_p
          delta := d_iter/d_crit;
          f.rd  := EoS.f_rd(delta=delta, tau=tau);
          RES_p := d_iter*T*R*(1+delta*f.rd) - p;
        end if;
    end while;
    // Modelica.Utilities.Streams.print("setState_pTX total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setState_pTX did not converge, input was p=" + String(p) + " and T=" + String(T));

    state.p := p;
    state.T := T;
    state.d := d_iter;
    f.i  := EoS.f_i(delta=delta, tau=tau);
    f.it := EoS.f_it(delta=delta, tau=tau);
    f.r  := EoS.f_r(delta=delta, tau=tau);
    f.rt := EoS.f_rt(delta=delta, tau=tau);
    state.h :=   T*R*(1 + tau*(f.it + f.rt) + delta*f.rd);
    state.u :=   T*R*(tau*(f.it+f.rt));
    state.s :=     R*(tau*(f.it+f.rt) - (f.i+f.r));
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
    EoS.HelmholtzDerivs f;

    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    SaturationProperties sat;
    MassFraction x "vapour quality";

    Density d_min;
    Density d_max;
    Density d_iter;
    Temperature T_min;
    Temperature T_max;
    Temperature T_iter;
    Real RES_p;
    Real RES_h;
    Real dpdd;
    Real dpdT;
    Real dhdd;
    Real dhdT;
    Real det "determinant of Jacobi matrix";
    Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
    Real tolerance=1e-9
    "tolerance for sum of relative RES_p and relative RES_h";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(p >= p_trip, "setState_phX_error: pressure is lower than triple point pressure");
      assert(p <= p_crit, "setState_phX_error: pressure is higher than critical pressure");
      sat := setSat_p(p=p);
      assert(h >= sat.liq.h, "setState_phX_error: enthalpy is lower than saturated liquid enthalpy: this is single phase liquid");
      assert(h <= sat.vap.h, "setState_phX_error: enthalpy is higher than saturated vapor enthalpy: this is single phase vapor");
    else
      if ((p < p_crit) and (p >= p_trip)) then
        // two-phase possible, do simple check first
        sat.Tsat := Ancillary.saturationTemperature_p(p=p);
        tau := T_crit/sat.Tsat;
        sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
        delta := sat.liq.d/d_crit;
        f.it  := EoS.f_it(delta=delta, tau=tau);
        f.rt  := EoS.f_rt(delta=delta, tau=tau);
        f.rd  := EoS.f_rd(delta=delta, tau=tau);
        sat.liq.h := sat.Tsat*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

        sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
        delta := sat.vap.d/d_crit;
        f.it  := EoS.f_it(delta=delta, tau=tau);
        f.rt  := EoS.f_rt(delta=delta, tau=tau);
        f.rd  := EoS.f_rd(delta=delta, tau=tau);
        sat.vap.h := sat.Tsat*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

        if ((h > sat.liq.h - abs(0.02*sat.liq.h)) and (h < sat.vap.h + abs(0.02*sat.vap.h))) then
          // Modelica.Utilities.Streams.print("two-phase state or close to it, get saturation properties from EoS", "printlog.txt");
          sat := setSat_p(p=p);
        end if;

        if (h < sat.liq.h) then
          // Modelica.Utilities.Streams.print("single phase liquid", "printlog.txt");
          state.phase := 1;
          d_min := sat.liq.d;
          d_max  := 1.1*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
          d_iter := sat.liq.d;
          T_min := fluidLimits.TMIN;
          T_max := sat.Tsat;
          T_iter:= sat.Tsat;
        elseif (h > sat.vap.h) then
          // Modelica.Utilities.Streams.print("single phase vapor", "printlog.txt");
          state.phase := 1;
          d_min := fluidLimits.DMIN;
          d_max := sat.liq.d;
          d_iter := sat.liq.d;
          T_min := sat.Tsat;
          T_max := fluidLimits.TMAX;
          T_iter:= sat.Tsat;
        else
          // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
          state.phase := 2;
        end if;

      else
        // Modelica.Utilities.Streams.print("p>=p_crit or p<p_trip, only single phase possible", "printlog.txt");
        state.phase := 1;
        d_min := fluidLimits.DMIN;
        d_max  := 1.1*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
        d_iter := d_crit;
        T_min := fluidLimits.TMIN;
        T_max := fluidLimits.TMAX;
        T_iter:= T_crit;
      end if;
    end if;

    // phase and region determination finished !

    if (state.phase == 2) then
      // force two-phase, SaturationProperties are already known
      state.p := p;
      state.h := h;
      state.T := sat.Tsat;
      x := (h - sat.liq.h)/(sat.vap.h - sat.liq.h);
      state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
    else
      // force single-phase
      f := EoS.setHelmholtzDerivs(d=d_iter, T=T_iter, phase=1);
      RES_p := d_iter*T_iter*f.R*(1+f.delta*f.rd) - p;
      RES_h := T_iter*f.R*((1+f.delta*f.rd)+f.tau*(f.it+f.rt)) - h;
      //RES_p := (1+delta*f.rd) - p/(d_iter*R*T_iter);
      //RES_h := ((1+delta*f.rd)+tau*(f.it+f.rt)) - h/(R*T_iter);

      while (((abs(RES_p/p) + abs(RES_h/h)) > tolerance) and (iter<iter_max)) loop
        iter := iter+1;

        // calculate gradients with respect to density and temperature
        dpdd := T_iter*f.R*(1+2*f.delta*f.rd+f.delta^2*f.rdd);
        dpdT := d_iter*f.R*(1+f.delta*f.rd-f.delta*f.tau*f.rtd);
        dhdd := T_iter*f.R/d_iter*(f.tau*f.delta*f.rtd+f.delta*f.rd+f.delta^2*f.rdd);
        dhdT := f.R*(-f.tau^2*(f.itt+f.rtt)+(1+f.delta*f.rd-f.delta*f.tau*f.rtd));

        /* // calculate dimensionless gradients of RES_p and RES_h with respect to delta and tau
      dpdd := (delta*f.rdd + 1*f.rd);
      dpdT := (delta*f.rtd);
      dhdd := (delta*f.rdd + 1*f.rd) + tau*(0+f.rtd);
      dhdT := (delta*f.rtd) + 1*(f.it+f.rt) +tau*(f.itt+f.rtt);*/

        // calculate determinant of Jacobi matrix
        det := dpdd*dhdT-dpdT*dhdd;

        /* // print for debugging
      Modelica.Utilities.Streams.print(" ", "printlog.txt");
      Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      Modelica.Utilities.Streams.print("d_iter=" + String(d_iter) + " and T_iter=" + String(T_iter), "printlog.txt");
      Modelica.Utilities.Streams.print("delta=" + String(delta) + " and tau=" + String(tau), "printlog.txt");
      Modelica.Utilities.Streams.print("RES_p=" + String(RES_p) + " and RES_h=" + String(RES_h), "printlog.txt");
      Modelica.Utilities.Streams.print("dpdd=" + String(dpdd) + " and dpdT=" + String(dpdT), "printlog.txt");
      Modelica.Utilities.Streams.print("dhdd=" + String(dhdd) + " and dhdT=" + String(dhdT), "printlog.txt");
      Modelica.Utilities.Streams.print("det(J)=" + String(det), "printlog.txt"); */

        // calculate better d_iter and T_iter
        d_iter := d_iter - gamma/det*(+dhdT*RES_p -dpdT*RES_h);
        T_iter := T_iter - gamma/det*(-dhdd*RES_p +dpdd*RES_h);

        /* // calculate better delta and tau
      delta := delta - gamma/det*(+dhdT*RES_p -dpdT*RES_h);
      tau   := tau   - gamma/det*(-dhdd*RES_p +dpdd*RES_h);
      d_iter := delta*d_crit;
      T_iter := T_crit/tau;*/

        // check bounds
        d_iter := max(d_min,d_iter);
        d_iter := min(d_max,d_iter);
        T_iter := max(T_min,T_iter);
        T_iter := min(T_max,T_iter);

        // calculate new RES_p and RES_h
        f := EoS.setHelmholtzDerivs(d=d_iter, T=T_iter, phase=1);
        RES_p := d_iter*T_iter*f.R*(1+f.delta*f.rd) - p;
        RES_h := T_iter*f.R*((1+f.delta*f.rd)+f.tau*(f.it+f.rt)) - h;
        //RES_p := (1+delta*f.rd) - p/(d_iter*R*T_iter);
        //RES_h := ((1+delta*f.rd)+tau*(f.it+f.rt)) - h/(R*T_iter);
      end while;
      // Modelica.Utilities.Streams.print("setState_phX total iteration steps " + String(iter), "printlog.txt");
      assert(iter<iter_max, "setState_phX did not converge, input was p=" + String(p) + " and h=" + String(h));

      state.p := p;
      state.h := h;
      state.d := d_iter;
      state.T := T_iter;
      state.u := T_iter*f.R*(f.tau*(f.it+f.rt));
      state.s :=        f.R*(f.tau*(f.it+f.rt) - (f.i+f.r));
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
    EoS.HelmholtzDerivs f;

    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    SaturationProperties sat;
    MassFraction x "vapour quality";

    Density d_min;
    Density d_max;
    Density d_iter;
    Temperature T_min;
    Temperature T_max;
    Temperature T_iter;
    Real RES_p;
    Real RES_s;
    Real dpdd;
    Real dpdT;
    Real dsdd;
    Real dsdT;
    Real det "determinant of Jacobi matrix";
    Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
    Real tolerance=1e-9
    "tolerance for sum of relative RES_p and relative RES_s ";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(p >= p_trip, "setState_psX_error: pressure is lower than triple point pressure");
      assert(p <= p_crit, "setState_psX_error: pressure is higher than critical pressure");
      sat := setSat_p(p=p);
      assert(s >= sat.liq.s, "setState_psX_error: entropy is lower than saturated liquid entropy: this is single phase liquid");
      assert(s <= sat.vap.s, "setState_psX_error: entropy is higher than saturated vapor entropy: this is single phase vapor");
    else
      if ((p < p_crit) and (p >= p_trip)) then
        // two-phase possible, do simple check first
        sat.Tsat := Ancillary.saturationTemperature_p(p=p);
        tau := T_crit/sat.Tsat;
        sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
        delta := sat.liq.d/d_crit;
        f.i   := EoS.f_i(tau=tau, delta=delta);
        f.it  := EoS.f_it(tau=tau, delta=delta);
        f.r   := EoS.f_r(tau=tau, delta=delta);
        f.rt  := EoS.f_rt(tau=tau, delta=delta);
        sat.liq.s := R*(tau*(f.it + f.rt) - f.i - f.r);

        sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
        delta := sat.vap.d/d_crit;
        f.i   := EoS.f_i(tau=tau, delta=delta);
        f.it  := EoS.f_it(tau=tau, delta=delta);
        f.r   := EoS.f_r(tau=tau, delta=delta);
        f.rt  := EoS.f_rt(tau=tau, delta=delta);
        sat.vap.s := R*(tau*(f.it + f.rt) - f.i - f.r);

        if ((s > sat.liq.s - abs(0.05*sat.liq.s)) and (s < sat.vap.s + abs(0.05*sat.vap.s))) then
          // Modelica.Utilities.Streams.print("two-phase state or close to it, get saturation properties from EoS", "printlog.txt");
          sat := setSat_p(p=p);
        end if;

        if (s < sat.liq.s) then
          // Modelica.Utilities.Streams.print("single phase liquid", "printlog.txt");
          state.phase := 1;
          d_min := sat.liq.d;
          d_max  := 1.1*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
          d_iter := sat.liq.d;
          T_min := fluidLimits.TMIN;
          T_max := sat.Tsat;
          T_iter:= sat.Tsat;
        elseif (s > sat.vap.s) then
          // Modelica.Utilities.Streams.print("single phase vapor", "printlog.txt");
          state.phase := 1;
          d_min := fluidLimits.DMIN;
          d_max := sat.liq.d;
          d_iter := sat.liq.d;
          T_min := sat.Tsat;
          T_max := fluidLimits.TMAX;
          T_iter:= sat.Tsat;
        else
          // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
          state.phase := 2;
        end if;

      else
        // Modelica.Utilities.Streams.print("p>=p_crit or p<p_trip, only single phase possible", "printlog.txt");
        state.phase := 1;
        d_min := fluidLimits.DMIN;
        d_max  := 1.1*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
        d_iter := d_crit;
        T_min := fluidLimits.TMIN;
        T_max := fluidLimits.TMAX;
        T_iter:= T_crit;
      end if;
    end if;

    // phase and region determination finished !

    if (state.phase == 2) then
      // force two-phase
      state.p := p;
      state.s := s;
      state.T := sat.Tsat;
      x := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
      state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
      state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
    else
      // force single-phase
      f := EoS.setHelmholtzDerivs(d=d_iter, T=T_iter, phase=1);
      RES_p := d_iter*T_iter*f.R*(1+f.delta*f.rd) - p;
      RES_s := f.R*(f.tau*(f.it + f.rt) - f.i - f.r) - s;

      while (((abs(RES_p/p) + abs(RES_s/s)) > tolerance) and (iter<iter_max)) loop
        iter := iter+1;

        // calculate gradients with respect to density and temperature
        dpdd := T_iter*f.R*(1+2*f.delta*f.rd+f.delta^2*f.rdd);
        dpdT := d_iter*f.R*(1+f.delta*f.rd-f.delta*f.tau*f.rtd);
        dsdd := f.R/d_iter*(-(1+f.delta*f.rd-f.delta*f.tau*f.rtd));
        dsdT := f.R/T_iter*(-f.tau^2*(f.itt+f.rtt));

        // calculate determinant of Jacobi matrix
        det := dpdd*dsdT-dpdT*dsdd;

        /* // print for debugging
      Modelica.Utilities.Streams.print(" ", "printlog.txt");
      Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      Modelica.Utilities.Streams.print("d_iter=" + String(d_iter) + " and T_iter=" + String(T_iter), "printlog.txt");
      Modelica.Utilities.Streams.print("RES_p=" + String(RES_p) + " and RES_s=" + String(RES_s), "printlog.txt");
      Modelica.Utilities.Streams.print("dpdd=" + String(dpdd) + " and dpdT=" + String(dpdT), "printlog.txt");
      Modelica.Utilities.Streams.print("dsdd=" + String(dsdd) + " and dsdT=" + String(dsdT), "printlog.txt");
      Modelica.Utilities.Streams.print("det(J)=" + String(det), "printlog.txt"); */

        // calculate better d_iter and T_iter
        d_iter := d_iter - gamma/det*(+dsdT*RES_p -dpdT*RES_s);
        T_iter := T_iter - gamma/det*(-dsdd*RES_p +dpdd*RES_s);

        // check bounds
        d_iter := max(d_min,d_iter);
        d_iter := min(d_max,d_iter);
        T_iter := max(T_min,T_iter);
        T_iter := min(T_max,T_iter);

        // calculate new RES_p and RES_s
        f := EoS.setHelmholtzDerivs(d=d_iter, T=T_iter, phase=1);
        RES_p := d_iter*T_iter*f.R*(1+f.delta*f.rd) - p;
        RES_s := f.R*(f.tau*(f.it + f.rt) - f.i - f.r) - s;
      end while;
      // Modelica.Utilities.Streams.print("setState_phX total iteration steps " + String(iter), "printlog.txt");
      assert(iter<iter_max, "setState_psX did not converge, input was p=" + String(p) + " and s=" + String(s));

      state.p := p;
      state.s := s;
      state.d := d_iter;
      state.T := T_iter;
      state.h := state.T*R*(f.tau*(f.it+f.rt) + (1+f.delta*f.rd));
      state.u := state.T*R*(f.tau*(f.it+f.rt));
    end if;

  end setState_psX;


  function setState_ThX "Return thermodynamic state as function of (T, h)"
    extends Modelica.Icons.Function;
    input Temperature T "Temperature";
    input SpecificEnthalpy h "Enthalpy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output ThermodynamicState state "thermodynamic state record";
protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Real delta "reduced density";
    Real tau(unit="1")=T_crit/T "inverse reduced temperature";
    EoS.HelmholtzDerivs f;

    SaturationProperties sat;
    MassFraction x "vapour quality";
    Density dmin=fluidLimits.DMIN;
    Density dmax=fluidLimits.DMAX;
    Real tolerance=1e-9 "relative Tolerance for Density";

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(T >= T_trip, "setState_ThX_error: pressure is lower than triple point pressure");
      assert(T <= T_crit, "setState_ThX_error: pressure is higher than critical pressure");
      sat := setSat_T(T=T);
      assert(h >= sat.liq.h, "setState_ThX_error: enthalpy is lower than saturated liquid enthalpy: this is single phase liquid");
      assert(h <= sat.vap.h, "setState_ThX_error: enthalpy is higher than saturated vapor enthalpy: this is single phase vapor");
    else
      if ((T <= T_crit) and (T >= T_trip)) then
        // two-phase possible, do simple check first
        sat.Tsat := T;
        tau := T_crit/sat.Tsat;
        sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
        delta := sat.liq.d/d_crit;
        f.it  := EoS.f_it(tau=tau, delta=delta);
        f.rt  := EoS.f_rt(tau=tau, delta=delta);
        f.rd  := EoS.f_rd(tau=tau, delta=delta);
        sat.liq.h := sat.Tsat*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

        sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
        delta := sat.vap.d/d_crit;
        f.it  := EoS.f_it(tau=tau, delta=delta);
        f.rt  := EoS.f_rt(tau=tau, delta=delta);
        f.rd  := EoS.f_rd(tau=tau, delta=delta);
        sat.vap.h := sat.Tsat*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

        if ((h > sat.liq.h - abs(0.05*sat.liq.h)) and (h < sat.vap.h + abs(0.05*sat.vap.h))) then
          // two-phase state or close to it, get saturation properties from EoS
          sat := setSat_T(T=sat.Tsat);
        end if;

        if (h < sat.liq.h) then
          state.phase := 1; // single phase liquid
          dmin := sat.liq.d;
        elseif (h > sat.vap.h) then
          state.phase := 1; // single phase vapor
          dmax := sat.vap.d;
        else
          state.phase := 2; // two-phase, all properties can be calculated from sat record
        end if;

      else
        // T>T_crit or T<T_trip, only single phase possible, do not change dmin and dmax
        state.phase := 1;
      end if;
    end if;

    state.T := T;
    state.h := h;
    if (state.phase == 2) then
      // force two-phase, SaturationProperties are already known
      state.p := sat.psat;
      x := (h - sat.liq.h)/(sat.vap.h - sat.liq.h);
      state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
    else
      // force single-phase
      state.d := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function setState_ThX_RES(
              T=T,
              h=h,
              phase=1),
            u_min=0.98*dmin,
            u_max=1.02*dmax,
            tolerance=tolerance);

      tau := T_crit/state.T;
      delta := state.d/d_crit;

      f.i   := EoS.f_i(tau=tau, delta=delta);
      f.it  := EoS.f_it(tau=tau, delta=delta);
      f.r   := EoS.f_r(tau=tau, delta=delta);
      f.rt  := EoS.f_rt(tau=tau, delta=delta);
      f.rd  := EoS.f_rd(tau=tau, delta=delta);
      state.p := state.d*T*R*(1+delta*f.rd);
      state.u := state.T*R*(tau*(f.it+f.rt));
      state.s :=         R*(tau*(f.it+f.rt) - (f.i+f.r));
    end if;

  end setState_ThX;


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
    assert(state.T >= T_trip, "vapourQuality error: Temperature is lower than triple-point temperature");
    assert(state.T <= T_crit, "vapourQuality error: Temperature is higher than critical temperature");

    if state.d <= sat.vap.d then
      x := 1;
    elseif state.d >= sat.liq.d then
      x := 0;
    else
      x := (1/state.d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
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
      f:=EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=1);
      cp := f.R*(-f.tau^2*(f.itt + f.rtt)
                 + (1 + f.delta*f.rd - f.delta*f.tau*f.rtd)^2/(1 + 2*f.delta*f.rd + f.delta^2*f.rdd));
  //  alternatively, yields same result
  //  cp := specificEnthalpy_derT_d(state,f) + specificEnthalpy_derd_T(state,f)*density_derT_p(state,f);
    elseif (state.phase == 2) then
      assert(false, "specificHeatCapacityCp warning: property not defined in two-phase region");
      cp := Modelica.Constants.inf; // division by zero
    end if;

  end specificHeatCapacityCp;


  redeclare function extends specificHeatCapacityCv
  "returns the isochoric specific heat capcacity"
  //input state
  //output cv

    // single phase
protected
    EoS.HelmholtzDerivs f;

    // two-phase
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta=state.d/d_crit "reduced density";
    Real tau=T_crit/state.T "inverse reduced temperature";

    SaturationProperties sat;
    DerPressureByTemperature dpT;
    EoS.HelmholtzDerivs
                    f_liq;
    EoS.HelmholtzDerivs
                    f_vap;
    DerPressureByTemperature dpTd_liq;
    DerPressureByTemperature dpTd_vap;
    DerPressureByDensity dpdT_liq;
    DerPressureByDensity dpdT_vap;
    SpecificHeatCapacity cv_lim2liq
    "limiting cv when approaching liq from within 2phase";
    SpecificHeatCapacity cv_lim2vap
    "limiting cv when approaching vap from within 2phase";
    MassFraction x "vapour quality";

  algorithm
    if (state.phase == 1) then
      f:=EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=1);
      cv := R*(-tau^2*(f.itt + f.rtt));
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T);
      // assert(false, "specificHeatCapacityCv warning: using cv in two-phase region", level=AssertionLevel.warning);
      // two-phase definition as in Span(2000), eq. 3.79 + 3.80 + 3.86
      // Attention: wrong sign in eq. 3.80
      dpT := saturationPressure_derT(T=state.T, sat=sat);
      f_liq := EoS.setHelmholtzDerivs(T=state.T, d=sat.liq.d, phase=1);
      f_vap := EoS.setHelmholtzDerivs(T=state.T, d=sat.vap.d, phase=1);
      dpTd_liq := pressure_derT_d(state=sat.liq);
      dpTd_vap := pressure_derT_d(state=sat.vap);
      dpdT_liq := pressure_derd_T(state=sat.liq);
      dpdT_vap := pressure_derd_T(state=sat.vap);

      cv_lim2liq := R*(-tau^2*(f_liq.itt + f_liq.rtt)) + state.T/sat.liq.d^2 * (dpTd_liq-dpT)^2/dpdT_liq;
      cv_lim2vap := R*(-tau^2*(f_vap.itt + f_vap.rtt)) + state.T/sat.vap.d^2 * (dpTd_vap-dpT)^2/dpdT_vap;

      x := (1/state.d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
      cv := cv_lim2liq + x*(cv_lim2vap-cv_lim2liq);
      // cv := (x*cv_lim2vap+(1-x)*cv_lim2liq) - (sat.vap.u-sat.liq.u)/(sat.vap.d-sat.liq.d)*(x*density_derT+(1-x)/sat.liq.d);
    end if;

  end specificHeatCapacityCv;


  redeclare function extends velocityOfSound
  "returns the speed or velocity of sound"
  //input state and output a are inherited from PartialMedium
  //input HelmholtzDerivs is optional and will be used for single-phase only
    input EoS.HelmholtzDerivs f=EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);

  algorithm
    assert(state.phase <> 2, "velocityOfSound error: property not defined in two-phase region");
    a := sqrt(state.T*f.R*( 1 + 2*f.delta*f.rd+ f.delta^2*f.rdd
              - ((1 + f.delta*f.rd) - f.delta*f.tau*f.rtd)^2 /
              (f.tau^2*(f.itt + f.rtt))));
  end velocityOfSound;


  redeclare function extends isobaricExpansionCoefficient
  "returns 1/v*(dv/dT)@p=const"
  //input state
  //output beta

protected
    EoS.HelmholtzDerivs f;
    Types.DerPressureByTemperature dpTd;
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      f:=EoS.setHelmholtzDerivs( T=state.T, d=state.d, phase=1);
      // Attention: wrong in Span(2000) table 3.10
      // correct in Lemmon(2000)
      // 1/v*(dv/dT)@p = -1/d*(dd/dT)@p = +1/d * (dp/dT)@d / (dp/dT)@d
      dpTd := pressure_derT_d(state=state);
      dpdT := pressure_derd_T(state=state);
      beta := 1/state.d*dpTd/dpdT;
    elseif (state.phase == 2) then
      beta := Modelica.Constants.small; // zero
    end if;
  end isobaricExpansionCoefficient;


  redeclare function extends isothermalCompressibility
  "returns -1/v*(dv/dp)@T=const"
  //input state
  //output kappa are inherited from PartialMedium

protected
    EoS.HelmholtzDerivs f;
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      f:=EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=1);
      dpdT := pressure_derd_T(state=state);
      kappa := 1/(state.d*dpdT);
    elseif (state.phase == 2) then
      kappa := Modelica.Constants.inf; // divide by zero
    end if;
  end isothermalCompressibility;


  redeclare replaceable function extends thermalConductivity
  "Return thermal conductivity"
    // inherits input state and output lambda
    // depends on dynamicViscosity, specificHeatCapacityCp, specificHeatCapacityCv and dpdd=1/dddp

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Density d_red_residual=fluidConstants[1].molarMass/
        thermalConductivityCoefficients.reducingMolarVolume_residual;
    Real delta "reduced density";

    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Temperature T_red_0=thermalConductivityCoefficients.reducingTemperature_0;
    Temperature T_red_residual=thermalConductivityCoefficients.reducingTemperature_residual;
    Real tau "reduced temperature";

    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    // coeffs for dilute contribution
    Real[size(thermalConductivityCoefficients.lambda_0_coeffs, 1),2] A=
        thermalConductivityCoefficients.lambda_0_coeffs;

    // coeffs for residual contribution
    Real[size(thermalConductivityCoefficients.lambda_r_coeffs, 1),4] B=
        thermalConductivityCoefficients.lambda_r_coeffs;

    // coeffs for critical enhancement
    Real nu=thermalConductivityCoefficients.nu;
    Real gamma=thermalConductivityCoefficients.gamma;
    Real R0=thermalConductivityCoefficients.R0;
    Real z=thermalConductivityCoefficients.z;
    Real c=thermalConductivityCoefficients.c;
    Real xi_0=thermalConductivityCoefficients.xi_0;
    Real Gamma_0=thermalConductivityCoefficients.Gamma_0;
    Real q_D=1/thermalConductivityCoefficients.qd_inverse;
    Temperature T_ref=thermalConductivityCoefficients.T_ref;

    // interim variables for critical enhancement
    constant Real pi=Modelica.Constants.pi;
    constant Real k_b=Modelica.Constants.k;
    Real ddpT;
    Real ddpT_ref;
    Real chi;
    Real chi_ref;
    Real Delta_chi;
    Real xi;
    Real Omega_0;
    Real Omega;

    SpecificHeatCapacity Cp;
    SpecificHeatCapacity Cv;
    DynamicViscosity eta_b;

    Real lambda_red_0=thermalConductivityCoefficients.reducingThermalConductivity_0;
    Real lambda_red_residual=thermalConductivityCoefficients.reducingThermalConductivity_residual;
    ThermalConductivity lambda_0=0;
    ThermalConductivity lambda_r=0;
    ThermalConductivity lambda_c=0;
    constant Real milli=1e-3;

  algorithm
    assert(state.phase <> 2, "thermalConductivity error: property not defined in two-phase region");

    // dilute gas contribution
    tau := state.T/T_red_0;
    lambda_0 := sum(A[i, 1]*tau^A[i, 2] for i in 1:size(A, 1));
    lambda_0 := lambda_0*lambda_red_0;

    // residual contribution; RefProp uses the name background contribution
    tau := state.T/T_red_residual;
    delta := state.d/d_red_residual;
    lambda_r := sum((B[i, 1]*tau^B[i, 2])*(delta)^B[i, 3] for i in 1:size(B, 1));
    lambda_r := lambda_r*lambda_red_residual;

    // critical enhancement by the simplified crossover model by Olchowy and Sengers
    if ((state.T > T_ref) or (state.d < d_crit/100)) then
      lambda_c := 0; // far away from critical point
    else
      // use critical values from EoS to calculate chi, Omega and lambda_c
      // watch out: algorithm for chi and chi_ref are different (chi_ref is multiplied with T_ref/state.T)
      ddpT := density_derp_T(state=state);
      chi := p_crit/d_crit^2*state.d*ddpT;

      delta := state.d/d_crit;
      tau := T_crit/T_ref;
      ddpT_ref := 1/(R*T_ref*(1 + 2*delta*EoS.f_rd(delta=delta, tau=tau) + delta^2*EoS.f_rdd(delta=delta, tau=tau)));
      chi_ref := p_crit/d_crit^2*state.d*ddpT_ref*T_ref/state.T;

      Delta_chi := chi - chi_ref;

      if (Delta_chi < 0) then
        lambda_c := 0;
      else
        xi := xi_0*(Delta_chi/Gamma_0)^(nu/gamma);

        Cp := specificHeatCapacityCp(state=state);
        Cv := specificHeatCapacityCv(state=state);
        Omega := 2/pi*((Cp - Cv)/Cp*atan(q_D*xi) + Cv/Cp*q_D*xi);
        Omega_0 := 2/pi*(1 - exp(-1/(1/(q_D*xi) + ((q_D*xi*d_crit/state.d)^2)/3)));

        eta_b := dynamicViscosity(state=state);
        lambda_c := (state.d*Cp*R0*k_b*state.T)/(6*pi*eta_b*xi)*(Omega - Omega_0);
        lambda_c := max(0, lambda_c);
      end if;
    end if;

    // RefPropresults are in mW/m·K but SI default is W/m·K
    lambda := milli*(lambda_0 + lambda_r + lambda_c);

    /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        d = " + String(state.d) + " and T = " + String(state.T));
  Modelica.Utilities.Streams.print(" lambda_0 = " + String(lambda_0));
  Modelica.Utilities.Streams.print(" lambda_r = " + String(lambda_r));
  Modelica.Utilities.Streams.print("   dpdd   = " + String(dpdd) + " and dpdd_ref = " + String(dpdd_ref));
  Modelica.Utilities.Streams.print("   chi    = " + String(chi) + "  and  chi_ref = " + String(chi_ref));
  Modelica.Utilities.Streams.print("Delta_chi = " + String(Delta_chi));
  Modelica.Utilities.Streams.print("       xi = " + String(xi));
  Modelica.Utilities.Streams.print("       Cp = " + String(Cp) + "  and  Cv = " + String(Cv));
  Modelica.Utilities.Streams.print("  Omega_0 = " + String(Omega_0));
  Modelica.Utilities.Streams.print("    Omega = " + String(Omega));
  Modelica.Utilities.Streams.print("    eta_b = " + String(eta_b));
  Modelica.Utilities.Streams.print(" lambda_c = " + String(lambda_c));
  Modelica.Utilities.Streams.print("  lambda  = " + String(lambda));
  Modelica.Utilities.Streams.print("===========================================");
  */

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
<dl>
<dt>Olchowy, G.A. and Sengers, J.V</dt>
<dd> <b>A simplified representation for the thermal conductivity of fluids in the critical region</b>.<br>
     International Journal of Thermophysics (1998) 10, 417-426.<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
  end thermalConductivity;


  redeclare replaceable function extends dynamicViscosity
  "Returns dynamic Viscosity"
    // inherits input state and output eta

protected
    DynamicViscosityModel dynamicViscosityModel=dynamicViscosityCoefficients.dynamicViscosityModel;
    CollisionIntegralModel collisionIntegralModel=dynamicViscosityCoefficients.collisionIntegralModel;
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Density d_red_residual=MM/dynamicViscosityCoefficients.reducingMolarVolume_residual;
    Real delta=0 "reduced density";
    Real delta_exp=0 "reduced density in exponential term";
    Real delta_0=0 "close packed density";
    Real dm=state.d/(1000*MM) "molar density in mol/l";     // 1 m3=1000 l
  //Real dm_crit=d_crit/(1000*MM) "molar density in mol/l"; // 1 m3=1000 l

    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Temperature T_red_0=dynamicViscosityCoefficients.reducingTemperature_0;
    Temperature T_red_residual=dynamicViscosityCoefficients.reducingTemperature_residual;
    Real T_star "reduced temperature";
    Real tau "reduced temperature";

    Real[size(dynamicViscosityCoefficients.a, 1),2] a=dynamicViscosityCoefficients.a;
    Real[size(dynamicViscosityCoefficients.b, 1),2] b=dynamicViscosityCoefficients.b;
    Real[size(dynamicViscosityCoefficients.c, 1),1] c=dynamicViscosityCoefficients.c;

    Real[size(dynamicViscosityCoefficients.g, 1),2] g=dynamicViscosityCoefficients.g;
    Real[size(dynamicViscosityCoefficients.e, 1),5] e=dynamicViscosityCoefficients.e;
    Real[size(dynamicViscosityCoefficients.nu_po, 1),5] nu_po=dynamicViscosityCoefficients.nu_po;
    Real[size(dynamicViscosityCoefficients.de_po, 1),5] de_po=dynamicViscosityCoefficients.de_po;
    // Real[size(dynamicViscosityCoefficients.nu_ex,1),5] nu_ex=dynamicViscosityCoefficients.nu_ex;
    // Real[size(dynamicViscosityCoefficients.de_ex,1),5] de_ex=dynamicViscosityCoefficients.de_ex;

    Real[size(dynamicViscosityCoefficients.CET, 1),2] CET=dynamicViscosityCoefficients.CET; // Chapman-Enskog-Term
    Real Omega=0 "reduced effective cross section / Omega collision integral";
    Real sigma=dynamicViscosityCoefficients.sigma;
    Real B_star=0 "reduced second viscosity virial coefficient";
    Real B=0 "second viscosity virial coefficient, l/mol";
    Real visci=0 "RefProp      visci temporary variable";
    Real xnum=0 "RefProp   numerator temporary variable";
    Real xden=0 "RefProp denominator temporary variable";
    Real G=0 "RefProp temporary variable";
    Real H=0 "RefProp temporary variable";
    Real F=0 "RefProp temporary variable";

    Real eta_red_0=dynamicViscosityCoefficients.reducingViscosity_0;
    Real eta_red_1=dynamicViscosityCoefficients.reducingViscosity_1;
    Real eta_red_residual=dynamicViscosityCoefficients.reducingViscosity_residual;
    DynamicViscosity eta_0=0 "zero density contribution";
    DynamicViscosity eta_1=0 "initial density contribution";
    DynamicViscosity eta_r=0 "residual viscosity";
    constant Real micro=1e-6;

  algorithm
    assert(state.phase <> 2, "dynamicViscosity error: property not defined in two-phase region");

    // collision integral
    if (collisionIntegralModel == CollisionIntegralModel.CI0) then
      T_star := (state.T/dynamicViscosityCoefficients.epsilon_kappa);
      Omega := 1.16145/T_star^0.14874 + 0.52487*exp(-0.77320*T_star) + 2.16178*exp(-2.43787*T_star);
    elseif (collisionIntegralModel == CollisionIntegralModel.CI1) then
      T_star := Modelica.Math.log(state.T/dynamicViscosityCoefficients.epsilon_kappa);
      Omega := exp(sum(a[i, 1]*(T_star)^a[i, 2] for i in 1:size(a, 1)));
    elseif (collisionIntegralModel == CollisionIntegralModel.CI2) then
      T_star := (dynamicViscosityCoefficients.epsilon_kappa/state.T)^(1/3);
      Omega := 1/(sum(a[i, 1]*(T_star)^(4-i) for i in 1:size(a, 1)));
    end if;

    // dilute gas (zero density) contribution
    // using the Chapman-Enskog-Term and the collision integral Omega
    if ((dynamicViscosityModel == DynamicViscosityModel.VS1)
    or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)) then
      tau := state.T/T_red_0;
      // first term is the Chapman-Enskog-Term
      eta_0 := CET[1, 1]*sqrt(tau)/(sigma^2*Omega);
      // possibly further empirical terms
      eta_0 := eta_0 + sum(CET[i, 1]*(tau)^CET[i, 2] for i in 2:size(CET, 1));
    elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
      eta_0 := CET[1, 1]*state.T^CET[1,2]/(sigma^2*Omega);
    elseif (dynamicViscosityModel == DynamicViscosityModel.VS4) then
    end if;
    eta_0 := eta_0*eta_red_0;

    // inital density contribution
    if ((dynamicViscosityModel == DynamicViscosityModel.VS1)
    or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)
    or  (dynamicViscosityModel == DynamicViscosityModel.VS4)) then
      // use the second viscosity virial coefficient B according to Rainwater and Friend
      T_star := (state.T/dynamicViscosityCoefficients.epsilon_kappa);
      B_star := sum(b[i, 1]*T_star^b[i, 2] for i in 1:size(b, 1));
      B := B_star*0.6022137*sigma^3;
      eta_1 := eta_0*B*dm;
    elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
      eta_1 := dm * (b[1,1] + b[2,1]*(b[3,1]-log(state.T/b[4,1]))^2);
    end if;
    eta_1 := eta_1*eta_red_1;

    // residual contribution
    if ((dynamicViscosityModel == DynamicViscosityModel.VS1)
    or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)) then
      // use the reduced close-packed density delta_0,
      // a simple polynominal, a rational polynominal and an exponential term
      tau := state.T/T_red_residual;
      delta := state.d/d_red_residual;
      if (abs(d_red_residual - 1) > 0.001) then
        delta_exp := state.d/d_crit;
      else
        delta_exp := delta;
      end if;

      if (dynamicViscosityModel == DynamicViscosityModel.VS1) then
        // generalized RefProp algorithm, be careful with coeffs: they may differ from article
        delta_0 := sum(g[i, 1]*tau^g[i, 2] for i in 1:size(g, 1));
      elseif (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative) then
        // alternative inverse form
        delta_0 := g[1, 1]/(1 + sum(g[i, 1]*tau^g[i, 2] for i in 2:size(g, 1)));
      end if;
      for i in 1:size(e, 1) loop
        visci := e[i, 1]*tau^e[i, 2]*delta^e[i, 3]*delta_0^e[i, 4]; // simple polynominal terms
        if (e[i, 5] > 0) then
          visci := visci*exp(-delta_exp^e[i, 5]);
        end if;
        eta_r := eta_r + visci;
      end for;

      for i in 1:size(nu_po, 1) loop
        // numerator of rational poly terms, RefProp algorithm
        xnum := xnum + (nu_po[i, 1]*tau^nu_po[i, 2]*delta^nu_po[i, 3]*delta_0^nu_po[i, 4]);
        if (nu_po[i, 5] > 0) then
          xnum := xnum*exp(-delta_exp^nu_po[i, 5]);
        end if;
      end for;
      for i in 1:size(de_po, 1) loop
        // denominator of rational poly terms, RefProp algorithm
        xden := xden + (de_po[i, 1]*tau^de_po[i, 2]*delta^de_po[i, 3]*delta_0^de_po[i, 4]);
        if (de_po[i, 5] > 0) then
          xden := xden*exp(-delta_exp^de_po[i, 5]);
        end if;
      end for;
      eta_r := eta_r + xnum/xden;
      // exponential terms not yet implemented!!

    elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
      G := c[1,1] + c[2,1]/state.T;
    //H := sqrt(dm)*(dm-dm_crit)/dm_crit;
      H := sqrt(dm)*(dm- c[8,1])/c[8,1];
      F := G + (c[3,1] + c[4,1]*state.T^(-3/2))*dm^0.1
             + (c[5,1] + c[6,1]/state.T + c[7,1]/state.T^2)*H;
      eta_r :=exp(F) - exp(G);

    elseif (dynamicViscosityModel == DynamicViscosityModel.VS4) then
      // not yet implemented!!
    end if;
    eta_r := eta_r*eta_red_residual;

    // RefProp results are in µPa·s where µ means micro or 1E-6 but SI default is Pa·s
    eta := micro*(eta_0 + eta_1 + eta_r);

    /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        T = " + String(state.T));
  Modelica.Utilities.Streams.print("   T_star = " + String(T_star));
  Modelica.Utilities.Streams.print("      tau = " + String(tau));
  Modelica.Utilities.Streams.print("        d = " + String(state.d));
  Modelica.Utilities.Streams.print("       dm = " + String(dm));
  Modelica.Utilities.Streams.print("    delta = " + String(delta));
  Modelica.Utilities.Streams.print("delta_exp = " + String(delta_exp));
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("    Omega = " + String(Omega));
  Modelica.Utilities.Streams.print("    eta_0 = " + String(eta_0));
  Modelica.Utilities.Streams.print("   B_star = " + String(B_star));
  Modelica.Utilities.Streams.print("        B = " + String(B));
  Modelica.Utilities.Streams.print("    eta_1 = " + String(eta_1));
  Modelica.Utilities.Streams.print("  delta_0 = " + String(delta_0));
  Modelica.Utilities.Streams.print("     xnum = " + String(xnum) + " and xden = " + String(xden));
  Modelica.Utilities.Streams.print("    eta_r = " + String(eta_r));
  Modelica.Utilities.Streams.print(" eta_r_RP = " + String(eta_r+eta_1));
  Modelica.Utilities.Streams.print("      eta = " + String(eta));
  Modelica.Utilities.Streams.print("===========================================");
  */

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


  redeclare function extends specificEnthalpy_dT
  // input, output and algorithm are inherited from PartialTwoPhaseMedium
  annotation (
    derivative=specificEnthalpy_dT_der);
  end specificEnthalpy_dT;


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


  function density_derT_h "returns density derivative (dd/dT)@h=const"
    input ThermodynamicState state "thermodynamic state record";
    output DerDensityByTemperature ddTh "Density derivative w.r.t. temperature";

protected
    SaturationProperties sat;
    Types.DerEnthalpyByTemperature dhTd;
    Types.DerEnthalpyByDensity dhdT;

  algorithm
    if (state.phase == 1) then
      dhTd := specificEnthalpy_derT_d(state=state);
      dhdT := specificEnthalpy_derd_T(state=state);
      ddTh := -dhTd/dhdT;
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T);
      ddTh := 50000000000000000;
    end if;
  end density_derT_h;


  redeclare function extends density_derp_T
  "returns density derivative (dd/dp)@T=const"
  //input state and output ddpT are inherited

protected
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      dpdT := pressure_derd_T(state=state);
      ddpT := 1.0/dpdT;
    elseif (state.phase == 2) then
      ddpT := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derp_T;


  redeclare function extends density_derT_p
  "returns density derivative (dd/dT)@p=const"
  //input state and output ddTp are inherited

protected
    Types.DerPressureByTemperature dpTd;
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      dpdT := pressure_derd_T(state=state);
      dpTd := pressure_derT_d(state=state);
      ddTp := -dpTd/dpdT;
    elseif (state.phase == 2) then
      ddTp := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derT_p;


  redeclare function extends density_derp_h
  "returns density derivative (dd/dp)@h=const"
  //input state
  //output ddph

protected
    SaturationProperties sat;
    Types.DerPressureByDensity dpdT;
    Types.DerPressureByTemperature dpTd;
          DerDensityByTemperature ddTh;

  algorithm
    if (state.phase == 1) then
      dpdT := pressure_derd_T(state=state);
      dpTd := pressure_derT_d(state=state);
      ddTh := density_derT_h(state=state);
      ddph := 1.0/(dpdT + dpTd/ddTh);
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T);
      dpTd := pressure_derT_d(state=state);
      ddph := state.d^2/state.T*specificHeatCapacityCv(state=state)/dpTd^2
              + state.d/state.T/dpTd;
    end if;
  end density_derp_h;


  redeclare function extends density_derh_p
  "returns density derivative (dd/dh)@p=const"
  //input state
  //output ddhp

protected
    SaturationProperties sat;
    Types.DerEnthalpyByDensity dhdT;
    Types.DerEnthalpyByTemperature dhTd;
          DerDensityByTemperature ddTp;

  algorithm
    if (state.phase == 1) then
      dhdT := specificEnthalpy_derd_T(state=state);
      dhTd := specificEnthalpy_derT_d(state=state);
      ddTp := density_derT_p(state=state);
      ddhp := 1.0/(dhdT + dhTd/ddTp);
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T);
      // dvhp = (v"-v')/(h"-h')
      // ddhp = -d^2 * dvhp
      ddhp := -state.d^2*(1/sat.liq.d-1/sat.vap.d)/(sat.liq.h-sat.vap.h);
    end if;
  end density_derh_p;


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

  input SaturationProperties sat=setSat_p(p=p) "optional input sat";
  // speeds up computation, if sat state is already known

  algorithm
    // inverse of (dp/dT)@sat
    // dTp := 1.0/saturationPressure_derT(T=sat.Tsat,sat=sat);
    // Clausius-Clapeyron, yields same result
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
    DerDensityByPressure ddpT = density_derp_T(state=sat.liq);
    DerDensityByTemperature ddTp = density_derT_p(state=sat.liq);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    ddldp := ddpT + ddTp*dTp;
  end dBubbleDensity_dPressure;


  redeclare function extends dDewDensity_dPressure
  "Return dew point density derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output ddvdp

protected
    DerDensityByPressure ddpT = density_derp_T(state=sat.vap);
    DerDensityByTemperature ddTp = density_derT_p(state=sat.vap);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    ddvdp := ddpT + ddTp*dTp;
  end dDewDensity_dPressure;


  redeclare function extends dBubbleEnthalpy_dPressure
  "Return bubble point enthalpy derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dhldp

protected
    DerEnthalpyByPressure dhpT = isothermalThrottlingCoefficient(state=sat.liq);
    DerEnthalpyByTemperature dhTp = specificHeatCapacityCp(state=sat.liq);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    dhldp := dhpT + dhTp*dTp;
  end dBubbleEnthalpy_dPressure;


  redeclare function extends dDewEnthalpy_dPressure
  "Return dew point enthalpy derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dhvdp

protected
    DerEnthalpyByPressure dhpT = isothermalThrottlingCoefficient(state=sat.vap);
    DerEnthalpyByTemperature dhTp = specificHeatCapacityCp(state=sat.vap);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    dhvdp := dhpT + dhTp*dTp;
  end dDewEnthalpy_dPressure;

end PartialHelmholtzMedium;
