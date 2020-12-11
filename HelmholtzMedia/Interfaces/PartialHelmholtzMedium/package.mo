within HelmholtzMedia.Interfaces;
partial package PartialHelmholtzMedium
  extends HelmholtzMedia.Interfaces.Types;


  extends HelmholtzMedia.Interfaces.Choices;


  extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(
    onePhase=false,
    singleState=false,
    smoothModel=true,
    DipoleMoment(min=0, max=5),
    AbsolutePressure(min=Modelica.Constants.small, max=1e12),
    SpecificEntropy(min=-Modelica.Constants.inf, max=Modelica.Constants.inf),
    ThermoStates=Choices.IndependentVariables.ph);

  constant HelmholtzMedia.Interfaces.Types.FluidLimits fluidLimits;

  constant EoS.HelmholtzCoefficients helmholtzCoefficients;

  constant Transport.ThermalConductivityCoefficients thermalConductivityCoefficients;

  constant Transport.DynamicViscosityCoefficients dynamicViscosityCoefficients;

  constant Transport.SurfaceTensionCoefficients surfaceTensionCoefficients;

  constant Ancillary.AncillaryCoefficients ancillaryCoefficients;

  // constant IndependentVariables independentVariables=IndependentVariables.dTX;
  constant InputChoice inputChoice=InputChoice.ph
  "Default choice of input variables for property computations";


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

  redeclare model extends BaseProperties(
      p(stateSelect = if preferredMediumStates and
                         (componentInputChoice == InputChoice.ph or
                          componentInputChoice == InputChoice.pT or
                          componentInputChoice == InputChoice.ps) then
                              StateSelect.prefer else StateSelect.default),
      T(stateSelect = if preferredMediumStates and
                         (componentInputChoice == InputChoice.pT or
                         componentInputChoice == InputChoice.dT) then
                           StateSelect.prefer else StateSelect.default),
      h(stateSelect = if preferredMediumStates and
                         componentInputChoice == InputChoice.ph then
                           StateSelect.prefer else StateSelect.default),
      d(stateSelect = if preferredMediumStates and
                         componentInputChoice == InputChoice.dT then
                           StateSelect.prefer else StateSelect.default))
  "Base properties (p, d, T, h, u, s) of a medium"

    SpecificEntropy s;
    parameter InputChoice componentInputChoice=inputChoice
    "Choice of input variables for property computations";

  equation
    MM = fluidConstants[1].molarMass;
    R = Modelica.Constants.R/MM;

    // use functions to calculate properties
    if (componentInputChoice == InputChoice.ph) then
      // state = setState_ph(p=p, h=h);
      d = density_ph(p=p, h=h);
      T = temperature_ph(p=p, h=h);
      s = specificEntropy_ph(p=p, h=h);
    elseif (componentInputChoice == InputChoice.ps) then
      // state = setState_ps(p=p, s=s);
      d = density_ps(p=p, s=s);
      T = temperature_ps(p=p, s=s);
      h = specificEnthalpy_ps(p=p, s=s);
    elseif (componentInputChoice == InputChoice.pT) then
      // state = setState_pT(p=p, T=T);
      d = density_pT(p=p, T=T);
      h = specificEnthalpy_pT(p=p, T=T);
      s = specificEntropy_pT(p=p, T=T);
    elseif (componentInputChoice == InputChoice.dT) then
      // state = setState_dT(d=d, T=T);
      p = pressure_dT(d=d, T=T);
      h = specificEnthalpy_dT(d=d, T=T);
      s = specificEntropy_dT(d=d, T=T);
    else
      assert(false, "Invalid choice for basePropertiesInput");
    end if;

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
    state.phase = if ((p < fluidConstants[1].criticalPressure) and (p > fluidConstants[1].triplePointPressure) and (h > sat.liq.h) and (h < sat.vap.h)) then 2 else 1;

  end BaseProperties;

  redeclare function setSat_T
  "iterative calculation of saturation properties from EoS with Newton-Raphson algorithm"
    // does not extend, because base class already defines an algorithm
    input Temperature T;
    output SaturationProperties sat;

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real tau(unit="1")=T_crit/T "inverse reduced temperature";

    EoS.HelmholtzDerivs fl(T=T);
    EoS.HelmholtzDerivs fv(T=T);

    Real delta_liq(unit="1", min=0);
    Real delta_vap(unit="1", min=0);
    Real J_liq;
    Real J_liq_delta;
    Real J_vap;
    Real J_vap_delta;
    Real K_liq;
    Real K_liq_delta;
    Real K_vap;
    Real K_vap_delta;

    Real RES[2] "residual function vector";
    Real Jacobian[2,2] "Jacobian matrix";
    Real NS[2] "Newton step vector";

    constant Real lambda(min=0.1,max=1) = 1 "convergence speed, default=1";
    constant Real tolerance=1e-7 "tolerance for p and d, needs to be smaller than simulation tolerance";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    // Modelica.Utilities.Streams.print("setSat_T: T="+String(T),"printlog.txt");
    sat.Tsat := T;

  if ((T>=T_trip) and (T<T_crit)) then
    // calculate guess values for reduced density delta
    delta_liq := 1.02*Ancillary.bubbleDensity_T(T=T)/d_crit;
    delta_vap := 0.98*Ancillary.dewDensity_T(T=T)/d_crit;

    // set necessary parts of fl and fv
    fl.r   := EoS.f_r(tau=tau, delta=delta_liq);
    fv.r   := EoS.f_r(tau=tau, delta=delta_vap);
    fl.rd  := EoS.f_rd(tau=tau, delta=delta_liq);
    fv.rd  := EoS.f_rd(tau=tau, delta=delta_vap);

    // dimensionless pressure
    J_liq := delta_liq*(1 + delta_liq*fl.rd);
    J_vap := delta_vap*(1 + delta_vap*fv.rd);
    // dimensionless Gibbs energy
    K_liq := delta_liq*fl.rd + fl.r + log(delta_liq);
    K_vap := delta_vap*fv.rd + fv.r + log(delta_vap);
    // residual vector
    RES := {(J_vap-J_liq), (K_vap-K_liq)};

    while (iter<iter_max) and (iter<1 or abs(RES[1])>tolerance or abs(RES[2])>tolerance or abs(NS[1])>tolerance or abs(NS[2])>tolerance) loop
      iter := iter+1;

      // calculate gradients of J and K, set up Jacobian matrix, get Newton step
      fl.rdd := EoS.f_rdd(tau=tau, delta=delta_liq);
      fv.rdd := EoS.f_rdd(tau=tau, delta=delta_vap);
      J_liq_delta := 1 + 2*delta_liq*fl.rd + delta_liq*delta_liq*fl.rdd;
      J_vap_delta := 1 + 2*delta_vap*fv.rd + delta_vap*delta_vap*fv.rdd;
      K_liq_delta := 2*fl.rd + delta_liq*fl.rdd + 1/delta_liq;
      K_vap_delta := 2*fv.rd + delta_vap*fv.rdd + 1/delta_vap;
      Jacobian := [-J_liq_delta, J_vap_delta;
                   -K_liq_delta, K_vap_delta];
      NS := -Modelica.Math.Matrices.solve(Jacobian,RES);

      // calculate better values for reduced density delta
      delta_liq := delta_liq + lambda*NS[1];
      delta_vap := delta_vap + lambda*NS[2];

      // check bounds
      delta_liq := max(delta_liq, 1);
      delta_liq := min(delta_liq, Modelica.Constants.inf);
      delta_vap := max(delta_vap, Modelica.Constants.small);
      delta_vap := min(delta_vap, 1);

      // update fl and fv
      fl.r   := EoS.f_r(tau=tau, delta=delta_liq);
      fv.r   := EoS.f_r(tau=tau, delta=delta_vap);
      fl.rd  := EoS.f_rd(tau=tau, delta=delta_liq);
      fv.rd  := EoS.f_rd(tau=tau, delta=delta_vap);

      // dimensionless pressure
      J_liq := delta_liq*(1 + delta_liq*fl.rd);
      J_vap := delta_vap*(1 + delta_vap*fv.rd);
      // dimensionless Gibbs energy
      K_liq := delta_liq*fl.rd + fl.r + log(delta_liq);
      K_vap := delta_vap*fv.rd + fv.r + log(delta_vap);
      // residual vector
      RES := {(J_vap-J_liq), (K_vap-K_liq)};

    end while;
    // Modelica.Utilities.Streams.print("setSat_T total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setSat_T did not converge, input was T=" + String(T) +
                          "; the remaining residuals are RES_J=" + String(RES[1]) +
                          " and RES_K=" + String(RES[2]));

    sat.liq := setState_dTX(d=delta_liq*d_crit, T=T, phase=1);
    sat.vap := setState_dTX(d=delta_vap*d_crit, T=T, phase=1);
    sat.psat := (sat.liq.p+sat.vap.p)/2;

  elseif (T>=T_crit) then
    // assert(T <= T_crit, "setSat_T error: Temperature is higher than critical temperature", level=AssertionLevel.warning);
    // above critical temperature, no stable two-phase state exists
    // anyway, it is possible to extend the vapour-pressure curve into this region
    // this can happen when called from BaseProperties
    // one possibility is use the state where ds/dT=max or ds/dp=max or dcp/dT=max or dcp/dp=max
    // here, the critical isochore is returned
    sat.liq := setState_dTX(d=d_crit, T=T, phase=1);
    sat.vap := setState_dTX(d=d_crit, T=T, phase=1);
    sat.psat := (sat.liq.p+sat.vap.p)/2;
  else
    // assert(T >= T_trip, "setSat_T error: Temperature is lower than triple-point temperature", level=AssertionLevel.warning);
    // T<T_trip: this does not make sense: if T is below the triple temperature, the medium is solid, not fluid
    // anyway, during initialization (at time=0) T=0 may happen
    // density values are extrapolated linearly, fantasy values are returned
    sat.Tsat := max(T, Modelica.Constants.small);
    delta_liq := (T_trip/sat.Tsat)*Ancillary.bubbleDensity_T(T=T_trip)/d_crit;
    delta_vap := (sat.Tsat/T_trip)*Ancillary.dewDensity_T(T=T_trip)/d_crit;

    sat.liq := setState_dTX(d=delta_liq*d_crit, T=T, phase=1);
    sat.vap := setState_dTX(d=delta_vap*d_crit, T=T, phase=1);
    sat.psat := (sat.liq.p+sat.vap.p)/2;
  end if;

  annotation (Documentation(info="
<html>
This function iteratively determines the saturation state  for a given temperature
by varying the density of saturated liquid and saturated vapor
with a Newton-Raphson approach for simultaneous equations.

<dl>
<dt> Ryo Akasaka:</dt>
<dd> <b>A reliable and useful method to determine the saturation state from Helmholtz energy equations of state</b>.<br />
     Journal of Thermal Science and Technology 3 (3) , 442-451.<br />
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
    constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
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
    Real Jacobian[3,3] "Jacobian matrix";
    Real NS[3] "Newton step vector";

    constant Real lambda(min=0.1,max=1) = 1 "convergence speed, default=1";
    constant Real tolerance=1e-6 "tolerance for T and d, needs to be smaller than simulation tolerance";
    Integer iter = 0;
    constant Integer iter_max = 200;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("setSat_p: p=" + String(p), "printlog.txt");
    sat.psat := p;

  if ((p>p_trip) and (p<p_crit)) then
    // calculate start values, density should be outside of two-phase dome
    sat.Tsat := Ancillary.saturationTemperature_p(p=p);
    sat.liq.d := 1.02*Ancillary.bubbleDensity_T(T=sat.Tsat);
    sat.vap.d := 0.98*Ancillary.dewDensity_T(T=sat.Tsat);

    // calculate residuals: RES=calc-input
    fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
    fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
    RES := {EoS.p(fl)-p, EoS.p(fv)-p, EoS.g(fl)-EoS.g(fv)};

    while (iter<iter_max) and (iter<1 or abs(RES[1])>tolerance or abs(RES[2])>tolerance or abs(NS[1])>tolerance or abs(NS[2])>tolerance or abs(NS[3])>tolerance) loop
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
    end while;
    // if verbose then Modelica.Utilities.Streams.print("setSat_p total iteration steps " + String(iter), "printlog.txt"); end if;
    // Modelica.Utilities.Streams.print("setSat_p total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setSat_p did not converge, input was p=" + String(p) +
                          "; the remaining residuals are RES_pl=" + String(RES[1]) +
                          " and RES_pv=" + String(RES[2]) +
                          " and RES_g=" + String(RES[3]));

    // check bounds, more strict
    sat.liq.d := max(sat.liq.d, d_crit);
    sat.liq.d := min(sat.liq.d, dl_trip);
    sat.vap.d := max(sat.vap.d, dv_trip);
    sat.vap.d := min(sat.vap.d, d_crit);
    sat.Tsat  := max(sat.Tsat,  T_trip);
    sat.Tsat  := min(sat.Tsat,  T_crit);

    sat.liq := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);

  elseif (p>=p_crit) then
    // assert(p <= p_crit, "setSat_p error: pressure is higher than critical pressure", level=AssertionLevel.warning);
    // above critical pressure, no stable two-phase state exists
    // anyway, it is possible to extend the vapour-pressure curve into this region
    // this can happen when called from BaseProperties
    // one possibility is to use the state where ds/dT=max or ds/dp=max or dcp/dT=max or dcp/dp=max
    // here, the critical isochore is returned
    sat.psat  := p;
    sat.liq := setState_pd(p=p, d=d_crit, phase=1);
    sat.vap := sat.liq;
    sat.Tsat := sat.liq.T;
  elseif (p<=p_trip) then
    // at pressures below triple point pressure, only single-phase vapor is possible
    // just return triple-point values
    sat.psat  := p;
    sat.liq := setState_dT(d=dl_trip, T=T_trip, phase=1);
    sat.vap := setState_dT(d=dv_trip, T=T_trip, phase=1);
    sat.Tsat := T_trip;
  else
    assert(false, "setSat_p: this should not happen, check p");
  end if;

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

  redeclare function extends setSmoothState
  "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
    import Modelica.Media.Common.smoothStep;

  algorithm
    state := ThermodynamicState(
      p=smoothStep(x, state_a.p, state_b.p, x_small),
      h=smoothStep(x, state_a.h, state_b.h, x_small),
      d=density_ph(p=smoothStep(x, state_a.p, state_b.p, x_small),
                   h=smoothStep(x, state_a.h, state_b.h, x_small)),
      T=temperature_ph(p=smoothStep(x, state_a.p, state_b.p, x_small),
                       h=smoothStep(x, state_a.h, state_b.h, x_small)),
      s=specificEntropy_ph(p=smoothStep(x, state_a.p, state_b.p, x_small),
                           h=smoothStep(x, state_a.h, state_b.h, x_small)),
      u=state.h-state.p/state.d,
      phase=0);
  annotation (Inline=true);
  end setSmoothState;

  redeclare function extends setState_dTX
  "Return thermodynamic state as function of (d, T)"

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant Temperature T_trip=fluidConstants[1].triplePointTemperature;

    EoS.HelmholtzDerivs f;
    SaturationProperties sat;
    MassFraction x "vapour quality";

  algorithm
    state.phase := phase;

    if state.phase==0 then
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
    elseif state.phase==2 then
      assert(T <= T_crit, "setState_dTX_error: Temperature is higher than critical temperature");
      sat := setSat_T(T=T);
      assert(d >= sat.vap.d, "setState_dTX_error: density is lower than saturated vapor density: this is single phase vapor");
      assert(d <= sat.liq.d, "setState_dTX_error: density is higher than saturated liquid density: this is single phase liquid");
    end if;

    state.d := d;
    state.T := T;
    if state.phase==2 then
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
    constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    // constant Real Z_crit=1/(d_crit*R*T_crit/p_crit);

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
    // Modelica.Utilities.Streams.print("setState_pTX: p=" + String(p) + " and T=" + String(T), "printlog.txt");
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
        d_max  := 1.10*fluidLimits.DMAX;
        d_iter := 1.10*sat.liq.d;
        // d_iter := Ancillary.density_pT_Soave(T=T, p=p, psat=sat.psat);
        // d_iter := 1/(R*T_crit/p_crit*Z_crit^(1+min(0,1-T/T_crit)^(2/7))) "Rackett";
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

    assert((RES_min*RES_max<0), "setState_pTX: d_min and d_max did not bracket the root", level=AssertionLevel.warning);
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
        // Modelica.Utilities.Streams.print("d_iter out of bounds, fallback to bisection method, step=" + String(iter) + ", d_iter=" + String(d_iter) + ", input was p=" + String(p) + " and T=" + String(T), "printlog.txt");
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
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau "inverse reduced temperature";

    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    constant SpecificEnthalpy h_crit=fluidConstants[1].HCRIT0;

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
    Integer iter_max = 200;
    Real lambda(min=1e-3,max=1) = 1 "convergence speed, default=1";

    Boolean useLineSearch=helmholtzCoefficients.useLineSearch;
    Integer iterLineSearch = 0;
    Real RSS_ls;
    Real lambda_ls;
    constant Real lambda_min = 0.01 "minimum for convergence speed";
    Real lambda_temp = 1 "temporary variable for convergence speed";
    constant Real alpha(min=0,max=1)=1e-4;
    Real rhs1;
    Real rhs2;
    Real a;
    Real b;
    Real Discriminant;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("setState_phX: p=" + String(p) + " and h=" + String(h), "printlog.txt");
    state.phase := phase;

    if state.phase==2 then
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
          T_iter:= sat.Tsat*0.98;
          T_max := sat.Tsat*1.02;
        elseif (h > sat.vap.h) then
          // Modelica.Utilities.Streams.print("single phase vapor", "printlog.txt");
          state.phase := 1;

          d_min := fluidLimits.DMIN;
          d_iter:= sat.vap.d/10;
          d_max := sat.vap.d*1.02;

          T_min := sat.Tsat*0.98;
          T_iter:= sat.Tsat*1.02;
          T_max := fluidLimits.TMAX*1.10;
        else
          // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
          state.phase := 2;
        end if;

      else
        state.phase := 1;
        // Modelica.Utilities.Streams.print("p>=p_crit, only single-phase possible", "printlog.txt");
        if (h<=h_crit) then
          // Modelica.Utilities.Streams.print("h<=h_crit, single-phase super-critical liquid-like region", "printlog.txt");
          d_min := d_crit*0.99;
          d_iter:= fluidLimits.DMAX*0.85;
          d_max := fluidLimits.DMAX*1.1;

          T_min := fluidLimits.TMIN;
          T_max := 1.25*Ancillary.saturationTemperature_h_liq(h=h);
          T_iter:= (1*T_min + 5*T_max/1.25)/6;
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

    if state.phase==2 then
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

        // calculate Jacobian matrix, Newton Step vector, gradient and slope
        Jacobian := [EoS.dpdT(f), EoS.dpTd(f);
                     EoS.dhdT(f), EoS.dhTd(f)];
        NS := -Modelica.Math.Matrices.solve(Jacobian,RES) "-J^(-1)*F";
        grad := RES*Jacobian "F*J";
        slope := grad*NS;
        // Modelica.Utilities.Streams.print("  Jacobian=" + Modelica.Math.Matrices.toString(Jacobian) + "  NS=" + Modelica.Math.Vectors.toString(NS) + "  grad=" + Modelica.Math.Vectors.toString(grad) + "  slope=" + String(slope), "printlog.txt");
        assert(slope<0,"roundoff problem, input was p=" + String(p) + " and h=" + String(h));

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
        // instead, the linesearching / backtracking loop tries to find lambda such that RSS decreases
        while useLineSearch and (lambda>=lambda_min) and not (RSS<=(RSS_old+alpha*lambda*slope)) loop
          iterLineSearch := iterLineSearch+1;
          iter_max := iter_max+1;

          // decrease lambda
          if (iterLineSearch<2) then
            // first attempt
            lambda_temp := -slope/(2*(RSS-RSS_old-slope));
          else
            // subsequent attempts
            rhs1 := RSS   -RSS_old-lambda   *slope;
            rhs2 := RSS_ls-RSS_old-lambda_ls*slope;
            a := (           rhs1/lambda^2 -        rhs2/lambda_ls^2) / (lambda-lambda_ls);
            b := (-lambda_ls*rhs1/lambda^2 + lambda*rhs2/lambda_ls^2) / (lambda-lambda_ls);
            if (a==0) then
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
              // new lambda should be less or equal 0.5*previous lambda
              lambda_temp := if (lambda_temp>0.5*lambda) then 0.5*lambda else lambda_temp;
            end if;
          end if;
          // store values for subsequent linesearch attempt
          lambda_ls := lambda;
          RSS_ls := RSS;

          // new lambda should be greater or equal 0.1*previous lambda
          lambda := max({lambda_temp, 0.1*lambda});

          // calculate new d_iter and T_iter using partial Newton step (lambda<1)
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
          // Modelica.Utilities.Streams.print("    linesearch attempt " + String(iterLineSearch) + ": lambda= " + String(lambda) + " d_iter=" + String(d_iter) + " T_iter= " + String(T_iter) + " RSS= " + String(RSS), "printlog.txt");

        end while;
        // reset iterLineSearch and lambda
        iterLineSearch := 0;
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
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau "inverse reduced temperature";

    constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    constant SpecificEntropy s_crit=fluidConstants[1].SCRIT0;

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

    Real[2] RES "residual function vector";
    Real RSS "residual sum of squares (divided by 2)";
    Real RSS_old "residual sum of squares (divided by 2)";
    Real[2,2] Jacobian "Jacobian matrix";
    Real[2] NS "Newton step vector";
    Real grad[2] "gradient vector";
    Real slope;

    EoS.HelmholtzDerivs f_med;
    Real[2] RES_med "residual function vector";
    Real[2] xy;
    Real[2] xy_old;
    Real[2] xy_med;
    Real[2,2] Jacobian_med "Jacobian matrix";

    constant Real tolerance=1e-7 "tolerance for RSS";
    Integer iter = 0;
    Integer iter_max = 200;
    Real lambda(min=1e-3,max=1) = 1 "convergence speed, default=1";

    Boolean useLineSearch=helmholtzCoefficients.useLineSearch;
    Integer iterLineSearch = 0;
    Real RSS_ls;
    Real lambda_ls;
    constant Real lambda_min = 0.01 "minimum for convergence speed";
    Real lambda_temp = 1 "temporary variable for convergence speed";
    constant Real alpha(min=0,max=1)=1e-4;
    Real rhs1;
    Real rhs2;
    Real a;
    Real b;
    Real Discriminant;

    ThermodynamicState state_ref;

  algorithm
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    // Modelica.Utilities.Streams.print("setState_psX: p=" + String(p) + " and s=" + String(s), "printlog.txt");
    state.phase := phase;

    if state.phase==2 then
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

          d_min := sat.liq.d*0.98;
          d_max := fluidLimits.DMAX*1.10;
          d_iter := sat.liq.d*1.10;

          T_min := fluidLimits.TMIN*0.99;
          T_max := sat.Tsat*1.02;
          T_iter:= sat.Tsat*0.95;
        elseif (s > sat.vap.s) then
          // Modelica.Utilities.Streams.print("single phase vapor", "printlog.txt");
          state.phase := 1;

          d_min := fluidLimits.DMIN;
          d_max := sat.vap.d*1.02;
          d_iter:= sat.vap.d/10;

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
        if (s<=s_crit) then
          // Modelica.Utilities.Streams.print("s<=s_crit, single-phase super-critical liquid-like region", "printlog.txt");
          state_ref := setState_pT(p=p, T=T_crit);
          if (s<state_ref.s) then
            T_min := max(fluidLimits.TMIN*0.99, Ancillary.saturationTemperature_s_liq(s=s)*0.95);
            T_max := T_crit;
            T_iter := 0.95*T_crit;
          elseif (s>state_ref.s) then
            T_min := T_crit;
            T_max := fluidLimits.TMAX*1.1;
            T_iter := 1.05*T_crit;
          else
            T_min := 0.95*T_crit;
            T_max := 1.05*T_crit;
            T_iter := T_crit;
          end if;
          // T_iter := min(T_min*1.5, T_max);
          // T_iter := (T_min+T_max)/2;
          // T_iter := min(T_min*1.5, T_crit);
          // T_iter:= (T_min+T_crit)/2;
          // T_iter := T_crit;

          d_min := d_crit*0.99;
          d_max := fluidLimits.DMAX*1.1;
          d_iter := fluidLimits.DMAX*0.9;
          // d_iter := min({state_ref.d, fluidLimits.DMAX*0.7});
          // d_iter := d_crit*1.02;
        else
          // Modelica.Utilities.Streams.print("s>s_crit, single-phase super-critical vapour-like region", "printlog.txt");
          state_ref := setState_pd(p=p, d=d_crit);
          if (s>state_ref.s) then
            d_min := fluidLimits.DMIN*0.99;
            d_max := d_crit;
            d_iter := 0.95*d_crit;
          elseif (s<state_ref.s) then
            d_min := d_crit;
            d_max := fluidLimits.DMAX*1.1;
            d_iter := 1.05*d_crit;
          else
            d_min := 0.95*d_crit;
            d_max := 1.05*d_crit;
            d_iter := d_crit;
          end if;
          // d_iter:= d_crit;
          // d_iter := p/(R_s*T_crit);
          // d_iter:= (d_min+d_max)/2;

          T_min := T_crit*0.98;
          T_max := fluidLimits.TMAX*1.1;
          // T_iter := state_ref.T;
          // T_iter := min(T_crit*2, T_max);
          T_iter:= T_crit*1.3;
          // T_iter := (T_min+T_max)/2;
        end if;
      end if;
    end if;

    // phase and region determination finished !

    if state.phase==2 then
      // Modelica.Utilities.Streams.print("two-phase, SaturationProperties are already known", "printlog.txt");
      state.p := p;
      state.s := s;
      state.T := sat.Tsat;
      x := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
      state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
      state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
      state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    else
      // Modelica.Utilities.Streams.print("single-phase, use 2D Newton-Raphson, start with d_iter=" + String(d_iter) + " (d_min=" + String(d_min) + " and d_max=" + String(d_max) + ") and T_iter=" + String(T_iter)+ " (T_min=" + String(T_min) + " and T_max=" + String(T_max) + ")", "printlog.txt");
      f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
      RES := {EoS.p(f)-p, EoS.s(f)-s};
      RSS := RES*RES/2;

      while ((RSS>tolerance) and (iter<iter_max)) loop
        iter := iter+1;

        // calculate Jacobian matrix, Newton Step vector, gradient and slope
        Jacobian := [EoS.dpdT(f), EoS.dpTd(f);
                     EoS.dsdT(f), EoS.dsTd(f)];
        NS := -Modelica.Math.Matrices.solve(Jacobian,RES) "-J^(-1)*F";
        grad := RES*Jacobian "F*J";
        slope := grad*NS;
        // Modelica.Utilities.Streams.print("  Jacobian=" + Modelica.Math.Matrices.toString(Jacobian) + "  NS=" + Modelica.Math.Vectors.toString(NS) + "  grad=" + Modelica.Math.Vectors.toString(grad) + "  slope=" + String(slope), "printlog.txt");
        assert(slope<0,"roundoff problem, input was p=" + String(p) + " and s=" + String(s), level=AssertionLevel.warning);

        // store old d_iter, T_iter and RSS
        d_iter_old := d_iter;
        T_iter_old := T_iter;
        RSS_old := RSS;

        // calculate new d_iter and T_iter using full Newton step
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];

        /* // Babajee
      xy_old := {d_iter,T_iter};
      xy_med := xy_old + NS;
      xy_med[1] := max(d_iter, d_min);
      xy_med[1] := min(d_iter, d_max);
      xy_med[2] := max(T_iter, T_min);
      xy_med[2] := min(T_iter, T_max);
      f_med := EoS.setHelmholtzDerivsSecond(d=xy_med[1], T=xy_med[2], phase=1);
      RES_med := {EoS.p(f_med)-p, EoS.s(f_med)-s};
      Jacobian_med := [EoS.dpdT(f_med), EoS.dpTd(f_med);
                       EoS.dsdT(f_med), EoS.dsTd(f_med)];
      xy := xy_med - Modelica.Math.Matrices.inv(Jacobian)*RES_med;
      d_iter:=xy[1];
      T_iter:=xy[2]; */

        // check bounds
        d_iter := max(d_iter, d_min);
        d_iter := min(d_iter, d_max);
        T_iter := max(T_iter, T_min);
        T_iter := min(T_iter, T_max);
        /*if (d_iter>d_max) then
        NS := NS*(d_iter-d_iter_old)/(d_max-d_iter_old);
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];
      end if;
      if (d_iter<d_min) then
        NS := NS*(d_min-d_iter_old)/(d_iter-d_iter_old);
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];
      end if;
      if (T_iter>T_max) then
        NS := NS*(T_iter-T_iter_old)/(T_max-T_iter_old);
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];
      end if;
      if (T_iter<T_min) then
        NS := NS*(T_min-T_iter_old)/(T_iter-T_iter_old);
        d_iter := d_iter_old + NS[1];
        T_iter := T_iter_old + NS[2];
      end if;*/

        // calculate new residual vector and residual sum of squares
        f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
        RES := {EoS.p(f)-p, EoS.s(f)-s};
        RSS := RES*RES/2;
        // Modelica.Utilities.Streams.print("iter=" + String(iter) + " d_iter=" + String(d_iter) + " T_iter=" + String(T_iter) + " RES_p=" + String(RES[1]) + " RES_s=" + String(RES[2]) + " RSS=" + String(RSS), "printlog.txt");

        // if RSS is not decreasing fast enough, the full Newton step is not used
        // instead, the linesearching / backtracking loop tries to find lambda such that RSS decreases
        while useLineSearch and (lambda>=lambda_min) and not (RSS<=(RSS_old+alpha*lambda*slope)) loop
          iterLineSearch := iterLineSearch+1;
          iter_max := iter_max+1;

          // decrease lambda
          if (iterLineSearch<2) then
            // first attempt
            lambda_temp := -slope/(2*(RSS-RSS_old-slope));
          else
            // subsequent attempts
            rhs1 := RSS   -RSS_old-lambda   *slope;
            rhs2 := RSS_ls-RSS_old-lambda_ls*slope;
            a := (           rhs1/lambda^2 -        rhs2/lambda_ls^2) / (lambda-lambda_ls);
            b := (-lambda_ls*rhs1/lambda^2 + lambda*rhs2/lambda_ls^2) / (lambda-lambda_ls);
            if (a==0) then
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
              // new lambda should be less or equal 0.5*previous lambda
              lambda_temp := if (lambda_temp>0.5*lambda) then 0.5*lambda else lambda_temp;
            end if;
          end if;
          // store values for subsequent linesearch attempt
          lambda_ls := lambda;
          RSS_ls := RSS;

          // new lambda should be greater or equal 0.1*previous lambda
          lambda := max({lambda_temp, 0.1*lambda});

          // calculate new d_iter and T_iter using partial Newton step (lambda<1)
          d_iter := d_iter_old +lambda*NS[1];
          T_iter := T_iter_old +lambda*NS[2];

          // check bounds
          d_iter := max(d_iter, d_min);
          d_iter := min(d_iter, d_max);
          T_iter := max(T_iter, T_min);
          T_iter := min(T_iter, T_max);

          // calculate new residual vector and residual sum of squares
          f := EoS.setHelmholtzDerivsSecond(d=d_iter, T=T_iter, phase=1);
          RES := {EoS.p(f)-p, EoS.s(f)-s};
          RSS := RES*RES/2;
          // Modelica.Utilities.Streams.print("    linesearch attempt " + String(iterLineSearch) + ": lambda= " + String(lambda) + " d_iter=" + String(d_iter) + " T_iter= " + String(T_iter) + " RSS= " + String(RSS), "printlog.txt");

        end while;
        // reset iterLineSearch and lambda
        iterLineSearch := 0;
        lambda_temp := 1;
        lambda := 1;

      end while;
      // Modelica.Utilities.Streams.print("setState_ps total iteration steps " + String(iter), "printlog.txt");
      assert(iter<iter_max, "setState_psX did not converge, input was p=" + String(p) + " and s=" + String(s) + ", remaining residual is RES_p=" +String(RES[1]) + " and RES_s=" + String(RES[2]));

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
  annotation(Inline = true);
  end temperature;

  redeclare function extends density
  "returns density from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output d
  algorithm
    d := state.d;
  annotation(Inline = true);
  end density;

  redeclare function extends pressure
  "returns pressure from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output p
  algorithm
    p := state.p;
  annotation(Inline = true);
  end pressure;

  redeclare function extends specificInternalEnergy
  "returns specificEnergy from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output u
  algorithm
    u := state.u;
  annotation(Inline = true);
  end specificInternalEnergy;

  redeclare function extends specificEntropy
  "returns specificEntropy from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output h
  algorithm
    s := state.s;
  annotation(Inline = true);
  end specificEntropy;

  redeclare function extends specificEnthalpy
  "returns specificEnthalpy from given ThermodynamicState"
  // inherited from: PartialMedium
  // inherits input state and output h
  algorithm
    h := state.h;
  annotation(Inline = true);
  end specificEnthalpy;

  redeclare function vapourQuality "returns the vapour quality"

  input ThermodynamicState state;
  output MassFraction x;

  algorithm
  x := vapourQuality_sat(state=state, sat=setSat_T(state.T));
  annotation (Inline=true);
  end vapourQuality;

  redeclare function extends specificHeatCapacityCp
  "returns the isobaric specific heat capcacity"
  //input state
  //output cp

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      cp := EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
    elseif state.phase==2 then
      assert(false, "specificHeatCapacityCp warning: in the two-phase region cp is infinite", level=AssertionLevel.warning);
      cp := Modelica.Constants.inf;
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
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      cv := EoS.duTd(f);

    elseif state.phase==2 then
      sat := setSat_T(T=state.T);
      x := vapourQuality_sat(state=state, sat=sat);
      dpT := saturationPressure_derT_sat(sat=sat);

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
  // a^2 := (dp/dd)@s=const

protected
    EoS.HelmholtzDerivs f;

    SaturationProperties sat;
    MassFraction x "vapour quality";
    DerTemperatureByPressure dTp;
    EoS.HelmholtzDerivs fl;
    EoS.HelmholtzDerivs fv;

    DerEntropyByPressure dsp_liq;
    DerEntropyByPressure dsp_vap;
    DerDensityByPressure ddp_liq;
    DerDensityByPressure ddp_vap;
    DerVolumeByPressure dvp_liq;
    DerVolumeByPressure dvp_vap;
    DerVolumeByPressure dvps;
    DerFractionByPressure dxps;

  algorithm
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      a := sqrt(EoS.dpdT(f)-EoS.dpTd(f)*EoS.dsdT(f)/EoS.dsTd(f));
    elseif state.phase==2 then
      sat := setSat_T(T=state.T);
      x := vapourQuality_sat(state=state, sat=sat);
      dTp := saturationTemperature_derp_sat(sat=sat);

      fl := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.liq.d, phase=1);
      fv := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.vap.d, phase=1);

      dsp_liq := EoS.dsdT(fl)/EoS.dpdT(fl) + (EoS.dsTd(fl)-EoS.dsdT(fl)*EoS.dpTd(fl)/EoS.dpdT(fl))*dTp;
      dsp_vap := EoS.dsdT(fv)/EoS.dpdT(fv) + (EoS.dsTd(fv)-EoS.dsdT(fv)*EoS.dpTd(fv)/EoS.dpdT(fv))*dTp;
      ddp_liq := 1.0/EoS.dpdT(fl) - EoS.dpTd(fl)/EoS.dpdT(fl)*dTp;
      ddp_vap := 1.0/EoS.dpdT(fv) - EoS.dpTd(fv)/EoS.dpdT(fv)*dTp;
      dvp_liq := -1.0/sat.liq.d^2 * ddp_liq;
      dvp_vap := -1.0/sat.vap.d^2 * ddp_vap;
      dxps :=(x*dsp_vap + (1 - x)*dsp_liq)/(sat.liq.s - sat.vap.s);

      dvps := dvp_liq + dxps*(1.0/sat.vap.d-1.0/sat.liq.d) + x*(dvp_vap-dvp_liq);
      a := sqrt(-1.0/(state.d*state.d*dvps));
    end if;
  end velocityOfSound;

  redeclare function extends isobaricExpansionCoefficient
  "returns 1/v*(dv/dT)@p=const"
  //input state
  //output beta

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if state.phase==1 then
      f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      beta := 1.0/state.d*EoS.dpTd(f)/EoS.dpdT(f);
    elseif state.phase==2 then
      assert(false, "isobaricExpansionCoefficient warning: in the two-phase region beta is zero", level=AssertionLevel.warning);
      beta := Modelica.Constants.small;
    end if;
  end isobaricExpansionCoefficient;

  redeclare function extends isothermalCompressibility
  "returns -1/v*(dv/dp)@T=const"
  //input state and
  //output kappa are inherited from PartialMedium

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if state.phase==1 then
      f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      kappa := 1.0/state.d/EoS.dpdT(f);
    elseif state.phase==2 then
      assert(false, "isothermalCompressibility warning: in the two-phase region kappa is infinite", level=AssertionLevel.warning);
      kappa := Modelica.Constants.inf;
    end if;
  end isothermalCompressibility;

  redeclare function extends isentropicExponent "returns -v/p*(dp/dv)@s=const"
  // also known as isentropic expansion coefficient

  algorithm
    gamma := state.d/state.p*velocityOfSound(state)^2;
    annotation(Inline = true);
  end isentropicExponent;

  redeclare function extends isentropicEnthalpy "returns isentropic enthalpy"

  algorithm
    h_is := specificEnthalpy(setState_ps(p=p_downstream, s=specificEntropy(refState)));
    annotation(Inline = true);
  end isentropicEnthalpy;

  redeclare replaceable function extends dynamicViscosity
  "Returns dynamic Viscosity"
    // inherits input state and output eta
  algorithm
    assert(state.phase <> 2, "dynamicViscosity warning: property not defined in two-phase region", level=AssertionLevel.warning);
    eta := Transport.dynamicViscosity(state);
  end dynamicViscosity;

  redeclare replaceable function extends thermalConductivity
  "Return thermal conductivity"
    // inherits input state and output lambda
  algorithm
    assert(state.phase <> 2, "thermalConductivity warning: property not defined in two-phase region", level=AssertionLevel.warning);
    lambda := Transport.thermalConductivity(state);
  end thermalConductivity;

  redeclare replaceable function extends surfaceTension
  "Return surface tension sigma in the two phase region"
      // inherits input saturationProperties sat and output SurfaceTension sigma
protected
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
  algorithm
    assert(sat.Tsat >= T_trip, "surfaceTension error: Temperature is lower than triple-point temperature");
    assert(sat.Tsat <= T_crit, "surfaceTension error: Temperature is higher than critical temperature");
    sigma := Transport.surfaceTension(sat=sat);
  end surfaceTension;

  redeclare function extends bubbleEnthalpy
  "returns specificEnthalpy from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output hl
  algorithm
    hl := sat.liq.h;
  annotation(Inline = true);
  end bubbleEnthalpy;

  redeclare function extends dewEnthalpy
  "returns specificEnthalpy from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output hv
  algorithm
    hv := sat.vap.h;
  annotation(Inline = true);
  end dewEnthalpy;

  redeclare function extends dewEntropy
  "returns specificEntropy from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output sv
  algorithm
    sv := sat.vap.s;
  annotation(Inline = true);
  end dewEntropy;

  redeclare function extends bubbleEntropy
  "returns specificEntropy from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output sl
  algorithm
    sl := sat.liq.s;
  annotation(Inline = true);
  end bubbleEntropy;

  redeclare function extends dewDensity
  "returns density from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dv
  algorithm
    dv := sat.vap.d;
  annotation(Inline = true);
  end dewDensity;

  redeclare function extends bubbleDensity
  "returns density from given SaturationProperties"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dl
  algorithm
    dl := sat.liq.d;
  annotation(Inline = true);
  end bubbleDensity;

  redeclare function extends saturationTemperature

  algorithm
    T := saturationTemperature_sat(setSat_p(p=p));
  annotation(Inline = true);
  end saturationTemperature;

  redeclare function extends saturationTemperature_derp "returns (dT/dp)@sat"

protected
    constant MolarMass MM = fluidConstants[1].molarMass;
    constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  algorithm
    dTp := if p<p_crit then saturationTemperature_derp_sat(sat=setSat_p(p=p)) else 1.0/pressure_derT_d(state=setState_pd(p=p, d=d_crit));
  annotation(Inline = true);
  end saturationTemperature_derp;

  redeclare function saturationTemperature_derp_sat "returns (dT/dp)@sat"
  // does not extend, because base class already defines an algorithm
  input SaturationProperties sat;
  output DerTemperatureByPressure dTp;

  algorithm
    assert(sat.vap.s<>sat.liq.s, "vap and liq entropy identical", level=AssertionLevel.warning);
    // vap volume is larger than liq volume, numerator is positive
    // vap entropy is larger than liq entropy, denominator is positive
    // Clausius-Clapeyron equation
    dTp := (1.0/sat.vap.d-1.0/sat.liq.d)/max(Modelica.Constants.eps, (sat.vap.s-sat.liq.s));
  annotation(Inline = true);
  end saturationTemperature_derp_sat;

  redeclare function extends saturationPressure

  algorithm
    p := saturationPressure_sat(setSat_T(T=T));
  annotation(Inline = true);
  end saturationPressure;

  redeclare function extends dBubbleDensity_dPressure
  "Return bubble point density derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output ddldp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.liq.T);
    DerDensityByPressure ddpT = 1.0/EoS.dpdT(f);
    DerDensityByTemperature ddTp = -EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = saturationTemperature_derp(p=sat.psat);

  algorithm
    ddldp := ddpT + ddTp*dTp;
  annotation(Inline = true);
  end dBubbleDensity_dPressure;

  redeclare function extends dDewDensity_dPressure
  "Return dew point density derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output ddvdp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.vap.T);
    DerDensityByPressure ddpT = 1.0/EoS.dpdT(f);
    DerDensityByTemperature ddTp = -EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = saturationTemperature_derp(p=sat.psat);

  algorithm
    ddvdp := ddpT + ddTp*dTp;
  annotation(Inline = true);
  end dDewDensity_dPressure;

  redeclare function extends dBubbleEnthalpy_dPressure
  "Return bubble point enthalpy derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dhldp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.liq.T);
    DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
    DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = saturationTemperature_derp(p=sat.psat);

  algorithm
    dhldp := dhpT + dhTp*dTp;
  annotation(Inline = true);
  end dBubbleEnthalpy_dPressure;

  redeclare function extends dDewEnthalpy_dPressure
  "Return dew point enthalpy derivative"
  // inherited from: PartialTwoPhaseMedium
  // inherits input sat and output dhvdp

protected
    EoS.HelmholtzDerivs f = EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.vap.T);
    DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
    DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
    DerTemperatureByPressure dTp = saturationTemperature_derp(p=sat.psat);

  algorithm
    dhvdp := dhpT + dhTp*dTp;
  annotation(Inline = true);
  end dDewEnthalpy_dPressure;

  redeclare function density_ph "returns density for given p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Enthalpy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output Density d "density";

  algorithm
    d := density_ph_state(p=p, h=h, state=setState_ph(p=p, h=h, phase=phase));

  annotation (
    Inline=true,
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
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
      ddph := 1.0/(EoS.dpdT(f) - EoS.dpTd(f)*EoS.dhdT(f)/EoS.dhTd(f));
    elseif state.phase==2 then
      sat := setSat_T(T=state.T);
      x := vapourQuality_sat(state=state, sat=sat);
      dTp := saturationTemperature_derp_sat(sat=sat);

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
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
      ddhp := 1.0/(EoS.dhdT(f) - EoS.dhTd(f)*EoS.dpdT(f)/EoS.dpTd(f));
    elseif state.phase==2 then
      sat:=setSat_T(T=state.T);
      // dvhp = (v"-v')/(h"-h')
      // ddhp = -d^2 * dvhp
      ddhp := -state.d^2*(1/sat.liq.d-1/sat.vap.d)/(sat.liq.h-sat.vap.h);
    end if;
  end density_derh_p;

  redeclare function temperature_ph "returns temperature for given p and h"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Enthalpy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output Temperature T "Temperature";

  algorithm
    T := temperature_ph_state(p=p, h=h, state=setState_ph(p=p, h=h, phase=phase));

  annotation (
    Inline=true,
    inverse(h=specificEnthalpy_pT(p=p, T=T, phase=phase)));
  end temperature_ph;

  redeclare function density_pT "Return density from p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output Density d "Density";

  algorithm
    d := density_pT_state(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase));

  annotation (
    Inline=true,
    inverse(p=pressure_dT(d=d, T=T, phase=phase),
            T=temperature_pd(p=p, d=d, phase=phase)));
  end density_pT;

  redeclare function extends density_derp_T
  "returns density derivative (dd/dp)@T=const"
  //input state and output ddpT are inherited

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      ddpT := 1.0/EoS.dpdT(f);
    elseif state.phase==2 then
      ddpT := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derp_T;

  redeclare function extends density_derT_p
  "returns density derivative (dd/dT)@p=const"
  //input state and output ddTp are inherited

protected
    EoS.HelmholtzDerivs f;

  algorithm
    if state.phase==1 then
      f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
      ddTp := -EoS.dpTd(f)/EoS.dpdT(f);
    elseif state.phase==2 then
      ddTp := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derT_p;

  redeclare function specificEnthalpy_pT
  "returns specific enthalpy for given p and T"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output SpecificEnthalpy h "specific enthalpy";

  algorithm
    h := specificEnthalpy_pT_state(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase));

  annotation (
    Inline=true,
    inverse(T=temperature_ph(p=p, h=h, phase=phase)));
  end specificEnthalpy_pT;

  redeclare function pressure_dT
    extends Modelica.Icons.Function;
    input Density d "Density";
    input Temperature T "Temperature";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output AbsolutePressure p "pressure";

  algorithm
    p := pressure_dT_state(d=d, T=T, state=setState_dT(d=d, T=T, phase=phase));

  annotation (
    Inline=true,
    inverse(d=density_pT(p=p, T=T, phase=phase),
            T=temperature_pd(p=p, d=d, phase=phase)));
  end pressure_dT;

  redeclare function specificEnthalpy_dT
    extends Modelica.Icons.Function;
    input Density d "Density";
    input Temperature T "Temperature";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output SpecificEnthalpy h "Specific Enthalpy";

  algorithm
    h := specificEnthalpy_dT_state(d=d, T=T, state=setState_dT(d=d, T=T, phase=phase));

  annotation (
    Inline=true);
  end specificEnthalpy_dT;

  redeclare function temperature_ps "returns temperature for given p and d"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Entropy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output Temperature T "Temperature";

  algorithm
    T := temperature(setState_ps(p=p, s=s, phase=phase));

  annotation (
    inverse(p=pressure_Ts(T=T, s=s, phase=phase),
            s=specificEntropy_pT(p=p, T=T, phase=phase)));
  end temperature_ps;

  redeclare function specificEnthalpy_ps
  "returns specific enthalpy for a given p and s"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Entropy";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  //input ThermodynamicState state;
    output SpecificEnthalpy h "specific enthalpy";

  algorithm
    h := specificEnthalpy(setState_psX(p=p, s=s, phase=phase));

  annotation (
    inverse(s=specificEntropy_ph(p=p, h=h, phase=phase)));
  end specificEnthalpy_ps;
end PartialHelmholtzMedium;
