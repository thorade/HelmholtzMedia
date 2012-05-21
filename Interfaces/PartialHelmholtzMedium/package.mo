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
  constant IndependentVariables independentVariables=IndependentVariables.dTX;

  constant HelmholtzCoefficients helmholtzCoefficients;

  constant ThermalConductivityCoefficients thermalConductivityCoefficients;

  constant DynamicViscosityCoefficients dynamicViscosityCoefficients;

  constant SurfaceTensionCoefficients surfaceTensionCoefficients;

  constant AncillaryCoefficients ancillaryCoefficients;


  redeclare function setSat_T
  "iterative calculation of saturation properties from EoS with Newton-Raphson algorithm"
    input Temperature T;
    input FixedPhase phase=2;
    output SaturationProperties sat;

protected
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real tau=T_crit/T "inverse reduced temperature";
    Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
    Real R=Modelica.Constants.R/fluidConstants[1].molarMass
    "specific gas constant in J/kg.K";

    Real delta_liq;
    Real delta_vap;
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
    Real Delta;
    Real gamma=1
    "parameter to control convergence speed, default is 1, decrease if convergence fails";
    Real tolerance=1e-9 "Tolerance for sum of Delta_J and Delta_K";

  algorithm
  if (phase==2) then
    // print("sat_of_T_EoS: T="+String(T),"printlog.txt");
    assert(T >= T_trip, "setSat_T error: Temperature is lower than triple-point temperature");
    assert(T <= T_crit, "setSat_T error: Temperature is higher than critical temperature");

    // calculate guess values for reduced density delta
    delta_liq := bubbleDensity_T_ANC(T=T)/d_crit;
    delta_vap := dewDensity_T_ANC(T=T)/d_crit;

    // pressure difference liquid-vapor
    J_liq := delta_liq*(1 + delta_liq*f_rd(tau=tau, delta=delta_liq));
    J_vap := delta_vap*(1 + delta_vap*f_rd(tau=tau, delta=delta_vap));
    Delta_J := abs(J_liq - J_vap);

    // Gibbs energy difference liquid-vapor
    K_liq := delta_liq*f_rd(tau=tau, delta=delta_liq) + f_r(tau=tau, delta=delta_liq) + log(delta_liq);
    K_vap := delta_vap*f_rd(tau=tau, delta=delta_vap) + f_r(tau=tau, delta=delta_vap) + log(delta_vap);
    Delta_K := abs(K_liq - K_vap);

    while (abs(Delta_J) + abs(Delta_K) > tolerance) loop
      // print("delta_liq=" + String(delta_liq) + " and delta_vap=" + String(delta_vap), "printlog.txt");
      // print("Delta_J=" + String(Delta_J) + " and Delta_K=" + String(Delta_K), "printlog.txt");

      // calculate better values for reduced density delta using gradients
      J_liq_delta := 1 + 2*delta_liq*f_rd(tau=tau, delta=delta_liq) + delta_liq^2*f_rdd(tau=tau, delta=delta_liq);
      J_vap_delta := 1 + 2*delta_vap*f_rd(tau=tau, delta=delta_vap) + delta_vap^2*f_rdd(tau=tau, delta=delta_vap);

      K_liq_delta := 2*f_rd(tau=tau, delta=delta_liq) + delta_liq*f_rdd(tau=tau, delta=delta_liq) + 1/delta_liq;
      K_vap_delta := 2*f_rd(tau=tau, delta=delta_vap) + delta_vap*f_rdd(tau=tau, delta=delta_vap) + 1/delta_vap;

      Delta := J_vap_delta*K_liq_delta - J_liq_delta*K_vap_delta;

      delta_liq := delta_liq + gamma/Delta*((K_vap - K_liq)*J_vap_delta - (J_vap - J_liq)*K_vap_delta);
      delta_vap := delta_vap + gamma/Delta*((K_vap - K_liq)*J_liq_delta - (J_vap - J_liq)*K_liq_delta);

      // calculate new Delta_J and Delta_K
      J_liq := delta_liq*(1 + delta_liq*f_rd(tau=tau, delta=delta_liq));
      J_vap := delta_vap*(1 + delta_vap*f_rd(tau=tau, delta=delta_vap));
      Delta_J := abs(J_liq - J_vap);

      K_liq := delta_liq*f_rd(tau=tau, delta=delta_liq) + f_r(tau=tau, delta=delta_liq) + log(delta_liq);
      K_vap := delta_vap*f_rd(tau=tau, delta=delta_vap) + f_r(tau=tau, delta=delta_vap) + log(delta_vap);
      Delta_K := abs(K_liq - K_vap);
    end while;

    sat.Tsat := T;
    sat.liq := setState_dTX(
          d=delta_liq*d_crit,
          T=T,
          phase=1);
    sat.vap := setState_dTX(
          d=delta_vap*d_crit,
          T=T,
          phase=1);
    sat.psat := sat.liq.p;
  // elseif (phase==1) then
  // print("sat_of_T_EoS has been called from single-phase state","printlog.txt");
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
    input FixedPhase phase=2;
    input Temperature T_guess=saturationTemperature(p=p); // optional input with default value
    output SaturationProperties sat;

protected
    Temperature T_min=max(0.98*T_guess, fluidConstants[1].triplePointTemperature);
    Temperature T_max=min(1.02*T_guess, fluidConstants[1].criticalTemperature);
    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    Temperature T;
    Real tolerance=1e-9 "relative Tolerance for Temperature";

  algorithm
  if (phase==2) then
    assert(p >= p_trip, "setSat_p error: Pressure is lower than triple-point pressure");
    assert(p <= p_crit, "setSat_p error: Pressure is higher than critical pressure");

    T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function setSat_p_RES(p=p),
          u_min=T_min,
          u_max=T_max,
          tolerance=tolerance);

    sat := setSat_T(T=T);
  // else: do nothing
  end if;
  end setSat_p;


  redeclare function extends saturationPressure
  "ancillary function: calculate saturation pressure for a given Temperature"
    // inherits input T and output p

protected
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real tau=T_crit/T "inverse reduced temperature";
    Real T_theta=1-T/T_crit;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    Integer nPressureSaturation=size(ancillaryCoefficients.pressureSaturation, 1);
    Real[nPressureSaturation] n=ancillaryCoefficients.pressureSaturation[:, 1];
    Real[nPressureSaturation] theta=ancillaryCoefficients.pressureSaturation[:, 2];

  algorithm
    assert(T <= T_crit, "saturationPressure error: Temperature is higher than critical temperature");
    p := p_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation));

    // this is an ancillary forward function
    // the corresponding iterative backward function is saturationTemperature(p)
    annotation (inverse(T=saturationTemperature(p=p)), Documentation(info="<html>
      <p>
      This algorithm returns the saturation pressure as a function of Temperature: psat=psat(T).
      This type of vapor pressure equation was developed by W. Wagner.
      Because it cannot be solved for temperature analytically, 
      the inverse function Tsat=Tsat(p) has to find Tsat iteratively.
      </p>
      
      <dl>
      <dt>Wagner, W.</dt>
      <dd> <b>Eine mathematisch statistische Methode zum Aufstellen thermodynamischer Gleichungen - gezeigt am Beispiel der Dampfdruckkurve reiner fluider Stoffe.</b><br>
           Forschrittberichte der VDI Zeitschriften, Reihe 3, Nr. 39 (1974)
      </dd>
      </dl>
      </html>"));
  end saturationPressure;


  redeclare function extends saturationTemperature
  "ancillary iterative function: calculate saturation temperature for a given pressure by iterating the anciallry function"
    // inherits input p and output T

protected
    Temperature T_trip=fluidConstants[1].triplePointTemperature;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
    Real tolerance=1e-9 "relative Tolerance for Density";

  algorithm
    assert(p >= p_trip, "saturationTemperature error: Pressure is lower than triple-point pressure");
    assert(p <= p_crit, "saturationTemperature error: Pressure is higher than critical pressure");
    T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function saturationTemperature_RES(p=p),
          u_min=T_trip,
          u_max=T_crit,
          tolerance=tolerance);

    // this is an iterative backward function
    // the corresponding ancillary forward function is saturationPressure(T)
    annotation (inverse(p=saturationPressure(T=T)));
  end saturationTemperature;


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


  redeclare record extends ThermodynamicState
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
  "Base properties (p, d, T, h, u, R, MM and, if applicable, X and Xi) of a medium"
    SpecificEntropy s;

  equation
    MM = fluidConstants[1].molarMass;
    R = Modelica.Constants.R/MM;

  algorithm
    if (independentVariables==IndependentVariables.dTX) then
      // Modelica.Utilities.Streams.print("Calculate thermodynamic state from EoS using setState_dT");
      state :=setState_dTX(d=d, T=T);
      p :=state.p;
      h :=state.h;
      s :=state.s;
      u :=h - p/d;
    elseif (independentVariables==IndependentVariables.pT) then
      // Modelica.Utilities.Streams.print("Calculate thermodynamic state from EoS using setState_pT");
      state :=setState_pTX(p=p, T=T);
      d :=state.d;
      h :=state.h;
      s :=state.s;
      u :=h - p/d;
    elseif (independentVariables==IndependentVariables.ph) then
      // Modelica.Utilities.Streams.print("Calculate thermodynamic state from EoS using setState_ph");
      state :=setState_phX(p=p, h=h);
      //d := density_ph(p=p, h=h);
       d :=state.d;
      T :=state.T;
      s :=state.s;
      u :=h - p/d;
    end if;

  if (state.phase == 2) then
      sat :=setSat_T(T=state.T);
    else
      // sat.Tsat :=saturationTemperature(p=p);
      // sat.psat :=saturationPressure(T=T);
    end if;

  end BaseProperties;


  redeclare function extends setState_dTX
  "Return thermodynamic state as function of d, T and composition X or Xi"

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta(unit="1")=d/d_crit "reduced density";
    Real tau(unit="1")=T_crit/T "inverse reduced temperature";
    HelmholtzDerivs f;

    SaturationProperties sat;
    MassFraction Q "vapour quality";

  algorithm
    state.phase := phase;

    if (state.phase == 0) then
      //phase unknown, check phase first
      if (T < T_crit) then
        // two-phase possible, do simple density check
        if ((d > 1.02*bubbleDensity_T_ANC(T=T)) or (d < 0.98*dewDensity_T_ANC(T=T))) then
          state.phase := 1;
        else
          // two-phase state or close to it, get saturation properties from EoS
          sat := setSat_T(T=T);
          if ((d < sat.liq.d) and (d > sat.vap.d)) then
            state.phase := 2;
          else
            state.phase := 1;
          end if;
        end if;
      else
        // T>T_crit
        state.phase := 1;
      end if;
    elseif (state.phase == 2) then
      assert(T <= T_crit, "setState_dTX_error: Temperature is higher than critical temperature");
      sat := setSat_T(T=T);
      assert(d >= sat.vap.d, "setState_dTX_error: density is lower than saturated vapor density: this is single phase vapor");
      assert(d <= sat.liq.d, "setState_dTX_error: density is higher than saturated liquid density: this is single phase liquid");
    end if;

    state.T := T;
    state.d := d;
    if (state.phase == 2) then
      // force two-phase
      Q := (1/d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
      state.p := sat.psat;
      state.h := sat.liq.h + Q*(sat.vap.h - sat.liq.h);
      state.u := sat.liq.u + Q*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + Q*(sat.vap.s - sat.liq.s);
    else
      // force single-phase
      f.i   := f_i(tau=tau, delta=delta);
      f.it  := f_it(tau=tau, delta=delta);
      f.r   := f_r(tau=tau, delta=delta);
      f.rt  := f_rt(tau=tau, delta=delta);
      f.rd  := f_rd(tau=tau, delta=delta);
      state.p := d*T*R*(1+delta*f.rd);
      state.h :=   T*R*(tau*(f.it + f.rt) + (1+delta*f.rd));
      state.u :=   T*R*(tau*(f.it + f.rt));
      state.s :=     R*(tau*(f.it + f.rt) - f.i - f.r);
    end if;

  end setState_dTX;


  redeclare function extends setState_pTX
  "Return thermodynamic state as function of p, T and composition X or Xi"

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau=T_crit/T "inverse reduced temperature";
    HelmholtzDerivs f;

    Real tolerance=1e-9 "relative Tolerance for Density";

  algorithm
    assert(phase <> 2, "setState_pTX_error: pressure and temperature are not independent varaibles in two-phase state");
    state.phase := 1;

    state.p := p;
    state.T := T;
    state.d := density_pT(
          p=p,
          T=T,
          phase=1);
    delta := state.d/d_crit;

    f.i   := f_i(tau=tau, delta=delta);
    f.it  := f_it(tau=tau, delta=delta);
    f.r   := f_r(tau=tau, delta=delta);
    f.rt  := f_rt(tau=tau, delta=delta);
    f.rd  := f_rd(tau=tau, delta=delta);
    state.h :=   T*R*(1 + tau*(f.it + f.rt) + delta*f.rd);
    state.u :=   T*R*(tau*(f.it+f.rt));
    state.s :=     R*(tau*(f.it+f.rt) - (f.i+f.r));
  end setState_pTX;


  redeclare function extends setState_phX
  "Return thermodynamic state as function of p, h and composition X or Xi"

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau "inverse reduced temperature";
    HelmholtzDerivs f;

    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    SaturationProperties sat;
    MassFraction Q "vapour quality";
    Temperature Tmin=fluidLimits.TMIN;
    Temperature Tmax=fluidLimits.TMAX;
    Real tolerance=1e-9 "relative Tolerance for Density";

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(p >= p_trip, "setState_phX_error: pressure is lower than triple point pressure");
      assert(p <= p_crit, "setState_phX_error: pressure is higher than critical pressure");
      sat := setSat_p(p=p);
      assert(h >= sat.liq.h, "setState_phX_error: enthalpy is lower than saturated liquid enthalpy: this is single phase liquid");
      assert(h <= sat.vap.h, "setState_phX_error: enthalpy is higher than saturated vapor enthalpy: this is single phase vapor");
    else
      if ((p <= p_crit) and (p >= p_trip)) then
        // two-phase possible, do simple check first
        sat.Tsat := saturationTemperature(p=p);
        tau := T_crit/sat.Tsat;
        sat.liq.d := bubbleDensity_T_ANC(T=sat.Tsat);
        delta := sat.liq.d/d_crit;
        f.it  := f_it(tau=tau, delta=delta);
        f.rt  := f_rt(tau=tau, delta=delta);
        f.rd  := f_rd(tau=tau, delta=delta);
        sat.liq.h := sat.Tsat*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

        sat.vap.d := dewDensity_T_ANC(T=sat.Tsat);
        delta := sat.vap.d/d_crit;
        f.it  := f_it(tau=tau, delta=delta);
        f.rt  := f_rt(tau=tau, delta=delta);
        f.rd  := f_rd(tau=tau, delta=delta);
        sat.vap.h := sat.Tsat*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

        if ((h > sat.liq.h - abs(0.05*sat.liq.h)) and (h < sat.vap.h + abs(0.05*sat.vap.h))) then
          // two-phase state or close to it, get saturation properties from EoS, use Tsat as starting value
          sat := setSat_p(p=p, T_guess=sat.Tsat);
        end if;

        if (h < sat.liq.h) then
          state.phase := 1; // single phase liquid
          Tmax := sat.Tsat;
        elseif (h > sat.vap.h) then
          state.phase := 1; // single phase vapor
          Tmin := sat.Tsat;
        else
          state.phase := 2; // two-phase, all properties can be calculated from sat record
        end if;

      else
        // p>p_crit or p<p_trip, only single phase possible, do not change Tmin and Tmax
        state.phase := 1;
      end if;
    end if;

    state.p := p;
    state.h := h;
    if (state.phase == 2) then
      // force two-phase, SaturationProperties are already known
      state.T := sat.Tsat;
      Q := (h - sat.liq.h)/(sat.vap.h - sat.liq.h);
      state.d := 1/(1/sat.liq.d + Q*(1/sat.vap.d - 1/sat.liq.d));
      state.u := sat.liq.u + Q*(sat.vap.u - sat.liq.u);
      state.s := sat.liq.s + Q*(sat.vap.s - sat.liq.s);
    else
      // force single-phase
      state.T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function setState_phX_RES(
              p=p,
              h=h,
              phase=1),
            u_min=0.98*Tmin,
            u_max=1.02*Tmax,
            tolerance=tolerance);
      state.d := density_pT(
            p=p,
            T=state.T,
            phase=1);

      tau := T_crit/state.T;
      delta := state.d/d_crit;

      f.i   := f_i(tau=tau, delta=delta);
      f.it  := f_it(tau=tau, delta=delta);
      f.r   := f_r(tau=tau, delta=delta);
      f.rt  := f_rt(tau=tau, delta=delta);
      state.u := state.T*R*(tau*(f.it+f.rt));
      state.s :=         R*(tau*(f.it+f.rt) - (f.i+f.r));
    end if;

  end setState_phX;


  redeclare function extends setState_psX
  "Return thermodynamic state as function of p, s and composition X or Xi"

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau "inverse reduced temperature";
    HelmholtzDerivs f;

    AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
    AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    SaturationProperties sat;
    MassFraction Q "vapour quality";
    Temperature Tmin=fluidLimits.TMIN;
    Temperature Tmax=fluidLimits.TMAX;
    Real tolerance=1e-9 "relative Tolerance for Density";

  algorithm
    state.phase := phase;

    if (state.phase == 2) then
      assert(p >= p_trip, "setState_psX_error: pressure is lower than triple point pressure");
      assert(p <= p_crit, "setState_psX_error: pressure is higher than critical pressure");
      sat := setSat_p(p=p);
      assert(s >= sat.liq.s, "setState_psX_error: entropy is lower than saturated liquid entropy: this is single phase liquid");
      assert(s <= sat.vap.s, "setState_psX_error: entropy is higher than saturated vapor entropy: this is single phase vapor");
    else
      if ((p <= p_crit) and (p >= p_trip)) then
        // two-phase possible, do simple check first
        sat.Tsat := saturationTemperature(p=p);
        tau := T_crit/sat.Tsat;
        sat.liq.d := bubbleDensity_T_ANC(T=sat.Tsat);
        delta := sat.liq.d/d_crit;
        f.i   := f_i(tau=tau, delta=delta);
        f.it  := f_it(tau=tau, delta=delta);
        f.r   := f_r(tau=tau, delta=delta);
        f.rt  := f_rt(tau=tau, delta=delta);
        sat.liq.s := R*(tau*(f.it + f.rt) - f.i - f.r);

        sat.vap.d := dewDensity_T_ANC(T=sat.Tsat);
        delta := sat.vap.d/d_crit;
        f.i   := f_i(tau=tau, delta=delta);
        f.it  := f_it(tau=tau, delta=delta);
        f.r   := f_r(tau=tau, delta=delta);
        f.rt  := f_rt(tau=tau, delta=delta);
        sat.vap.s := R*(tau*(f.it + f.rt) - f.i - f.r);

        if ((s > sat.liq.s - abs(0.05*sat.liq.s)) and (s < sat.vap.s + abs(0.05*sat.vap.s))) then
          // two-phase state or close to it, get saturation properties from EoS, use Tsat as starting value
          sat := setSat_p(p=p, T_guess=sat.Tsat);
        end if;

        if (s < sat.liq.s) then
          state.phase := 1; // single phase liquid
          Tmax := sat.Tsat;
        elseif (s > sat.vap.s) then
          state.phase := 1; // single phase vapor
          Tmin := sat.Tsat;
        else
          state.phase := 2; // two-phase, all properties can be calculated from sat record
        end if;

      else
        // p>p_crit or p<p_trip, only single phase possible, do not change Tmin and Tmax
        state.phase := 1;
      end if;
    end if;

    state.p := p;
    state.s := s;
    if (state.phase == 2) then
      // force two-phase
      state.T := sat.Tsat;
      Q := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
      state.d := 1/(1/sat.liq.d + Q*(1/sat.vap.d - 1/sat.liq.d));
      state.h := sat.liq.h + Q*(sat.vap.h - sat.liq.h);
      state.u := sat.liq.u + Q*(sat.vap.u - sat.liq.u);
    else
      // force single-phase
      state.T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function setState_psX_RES(
              p=p,
              s=s,
              phase=1),
            u_min=0.98*Tmin,
            u_max=1.02*Tmax,
            tolerance=tolerance);
      state.d := density_pT(
            p=p,
            T=state.T,
            phase=1);

      tau := T_crit/state.T;
      delta := state.d/d_crit;

      f.it  := f_it(tau=tau, delta=delta);
      f.rt  := f_rt(tau=tau, delta=delta);
      f.rd  := f_rd(tau=tau, delta=delta);
      state.h := state.T*R*(tau*(f.it+f.rt) + (1+delta*f.rd));
      state.u := state.T*R*(tau*(f.it+f.rt));
    end if;

  end setState_psX;


  redeclare function density_pT
  "iteratively finds the density for a given p and T (works for single-phase only)"

    // this function will be called millions of times,
    // so it makes sense to have good starting values u_min and u_max
    // an estimate for rho can be calculated from Redlich-Kwong-Soave
    // then u_min=0.8*estimate and u_max=1.2*estimate

    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input FixedPhase phase=1 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Density d "Density";

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta=d/d_crit "reduced density";
    Real tau=T_crit/T "inverse reduced temperature";

    SaturationProperties sat;

    Density dmin=fluidLimits.DMIN;
    Density dmax=fluidLimits.DMAX;
    Real tolerance=1e-9 "relative Tolerance for Density";

  algorithm
    assert(phase <> 2, "density_pT error: in two-phase state pressure and temperature are not independent variables");

    if (T <= T_crit) then
      sat.psat := saturationPressure(T=T);

      if (p > 1.05*sat.psat) then
        sat.liq.d := bubbleDensity_T_ANC(T=T);
      elseif (p < 0.95*sat.psat) then
        sat.vap.d := dewDensity_T_ANC(T=T);
      else
        sat := setSat_T(T=T); // very close to two-phase: get saturation properties from EoS for consistency
      end if;

      if (p > sat.psat) then
        dmin := sat.liq.d; // single-phase liquid: d is between dliq and rho_max
      elseif (p < sat.psat) then
        dmax := sat.vap.d; //single-phase vapor: d is between 0 and dvap
      else
        assert(p <> sat.psat, "density_pT error: cannot calculate the density because in the two-phase region p and T are not independent");
      end if;

    end if;

    d := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
         function density_pT_RES(T=T, p=p),
          u_min=0.95*dmin,
          u_max=1.05*dmax,
          tolerance=tolerance);

    // this is an iterative backward function
    // pressure_dT is the corresponding forward function
     annotation (inverse(p=pressure_dT(d=d, T=T, phase=phase)));
  end density_pT;


  redeclare function specificEnthalpy_pT
  "iteratively finds the specific enthalpy for a given p and T"

    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output SpecificEnthalpy h "Specific Enthalpy";

protected
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta "reduced density";
    Real tau=T_crit/T "inverse reduced temperature";
    HelmholtzDerivs f;

    Density d;

  algorithm
    assert(phase <> 2, "specificEnthalpy_pT error: pressure and temperature are not independent variables in two-phase state");
    d := density_pT(p=p, T=T, phase=phase);
    delta := d/d_crit;
    f.it  := f_it(tau=tau, delta=delta);
    f.rt  := f_rt(tau=tau, delta=delta);
    f.rd  := f_rd(tau=tau, delta=delta);
    h :=   T*R*(1 + tau*(f.it + f.rt) + delta*f.rd);

    // this is an iterative backward function
    // the two inverse functions are T=temperature_ph and p=pressure_Th
    annotation (inverse(T=temperature_ph(p=p, h=h, phase=phase)));
  end specificEnthalpy_pT;


  redeclare function extends specificHeatCapacityCp
  "returns the isobaric specific heat capcacity"
  //input state
  //output cp

protected
    HelmholtzDerivs f;

  algorithm
    if (state.phase == 1) then
      f:=setHelmholtzDerivs(
          T=state.T,
          d=state.d,
          phase=1);
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
    HelmholtzDerivs f;

    // two-phase
    MolarMass MM = fluidConstants[1].molarMass;
    SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
    Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
    Temperature T_crit=fluidConstants[1].criticalTemperature;
    Real delta=state.d/d_crit "reduced density";
    Real tau=T_crit/state.T "inverse reduced temperature";

    SaturationProperties sat;
    DerPressureByTemperature dpT;
    HelmholtzDerivs f_liq;
    HelmholtzDerivs f_vap;
    DerPressureByTemperature dpTd_liq;
    DerPressureByTemperature dpTd_vap;
    DerPressureByDensity dpdT_liq;
    DerPressureByDensity dpdT_vap;
    SpecificHeatCapacity cv_lim2liq
    "limiting cv when approaching liq from within 2phase";
    SpecificHeatCapacity cv_lim2vap
    "limiting cv when approaching vap from within 2phase";
    MassFraction Q "vapour quality";

  algorithm
    if (state.phase == 1) then
      f:=setHelmholtzDerivs(
          T=state.T,
          d=state.d,
          phase=1);
      cv := R*(-tau^2*(f.itt + f.rtt));
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T, phase=state.phase);
      // assert(false, "specificHeatCapacityCv warning: using cv in two-phase region", level=AssertionLevel.warning);
      // two-phase definition as in Span(2000), eq. 3.79 + 3.80 + 3.86
      // Attention: wrong sign in eq. 3.80
      dpT := saturationPressure_derT(T=state.T, sat=sat);
      f_liq := setHelmholtzDerivs(T=state.T, d=sat.liq.d, phase=1);
      f_vap := setHelmholtzDerivs(T=state.T, d=sat.vap.d, phase=1);
      dpTd_liq := pressure_derT_d(state=sat.liq, f=f_liq);
      dpTd_vap := pressure_derT_d(state=sat.vap, f=f_vap);
      dpdT_liq := pressure_derd_T(state=sat.liq, f=f_liq);
      dpdT_vap := pressure_derd_T(state=sat.vap, f=f_vap);

      cv_lim2liq := R*(-tau^2*(f_liq.itt + f_liq.rtt)) + state.T/sat.liq.d^2 * (dpTd_liq-dpT)^2/dpdT_liq;
      cv_lim2vap := R*(-tau^2*(f_vap.itt + f_vap.rtt)) + state.T/sat.vap.d^2 * (dpTd_vap-dpT)^2/dpdT_vap;

      Q := (1/state.d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
      cv := cv_lim2liq + Q*(cv_lim2vap-cv_lim2liq);
      // cv := (Q*cv_lim2vap+(1-Q)*cv_lim2liq) - (sat.vap.u-sat.liq.u)/(sat.vap.d-sat.liq.d)*(Q*density_derT+(1-Q)/sat.liq.d);
    end if;

  end specificHeatCapacityCv;


  redeclare function extends velocityOfSound
  "returns the speed or velocity of sound"
  //input state and output a are inherited from PartialMedium
  //input HelmholtzDerivs is optional and will be used for single-phase only
    input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);

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
    HelmholtzDerivs f;
    Types.DerPressureByTemperature dpTd;
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      f:=setHelmholtzDerivs(
          T=state.T,
          d=state.d,
          phase=1);
      // Attention: wrong in Span(2000) table 3.10
      // correct in Lemmon(2000)
      // 1/v*(dv/dT)@p = -1/d*(dd/dT)@p = +1/d * (dp/dT)@d / (dp/dT)@d
      dpTd := pressure_derT_d(state=state, f=f);
      dpdT := pressure_derd_T(state=state, f=f);
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
    HelmholtzDerivs f;
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      f:=setHelmholtzDerivs(
          T=state.T,
          d=state.d,
          phase=1);
      dpdT := pressure_derd_T(state=state, f=f);
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
      ddpT_ref := 1/(R*T_ref*(1 + 2*delta*f_rd(delta=delta, tau=tau) + delta^2*f_rdd(delta=delta, tau=tau)));
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

    Real[1,2] CET=dynamicViscosityCoefficients.CET; // Chapman-Enskog-Term
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
      eta_0 := CET[1, 1]*sqrt(tau)/(sigma^2*Omega);
    elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
      eta_0 := CET[1, 1]*state.T^CET[1,2]/(sigma^2*Omega);
    elseif (dynamicViscosityModel == DynamicViscosityModel.VS4) then
    end if;
    eta_0 := eta_0*eta_red_0;

    // inital density contribution
    if ((dynamicViscosityModel == DynamicViscosityModel.VS1)
    or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)
    or  (dynamicViscosityModel == DynamicViscosityModel.VS4)) then
      // use the second viscosity virial coefficient B
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


  function density_derT_h "returns density derivative (dd/dT)@h=const"
    input ThermodynamicState state "thermodynamic state record";
  //input HelmholtzDerivs is optional and will be used for single-phase only
    input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
  //input sat is optional and will be used for two-phase only
    input SaturationProperties sat=setSat_T(T=state.T, phase=state.phase);
    output DerDensityByTemperature ddTh "Density derivative w.r.t. temperature";

protected
    Types.DerEnthalpyByTemperature dhTd;
    Types.DerEnthalpyByDensity dhdT;

  algorithm
    if (state.phase == 1) then
      dhTd := specificEnthalpy_derT_d(state=state, f=f);
      dhdT := specificEnthalpy_derd_T(state=state, f=f);
      ddTh := -dhTd/dhdT;
    elseif (state.phase == 2) then
      ddTh := 50000000000000000;
    end if;
  end density_derT_h;


  redeclare function extends density_derp_T
  "returns density derivative (dd/dp)@T=const"
  //input state and output ddpT are inherited
  //input HelmholtzDerivs is optional and will be used for single-phase only
    input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);

protected
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      dpdT := pressure_derd_T(state=state, f=f);
      ddpT := 1.0/dpdT;
    elseif (state.phase == 2) then
      ddpT := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derp_T;


  redeclare function extends density_derT_p
  "returns density derivative (dd/dT)@p=const"
  //input state and output ddTp are inherited
  //input HelmholtzDerivs is optional and will be used for single-phase only
    input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);

protected
    Types.DerPressureByTemperature dpTd;
    Types.DerPressureByDensity dpdT;

  algorithm
    if (state.phase == 1) then
      dpdT := pressure_derd_T(state=state, f=f);
      dpTd := pressure_derT_d(state=state, f=f);
      ddTp := -dpTd/dpdT;
    elseif (state.phase == 2) then
      ddTp := Modelica.Constants.inf; // divide by zero
    end if;
  end density_derT_p;


  redeclare function extends density_ph "returns density for given p and h"
  // inherited from: PartialMedium
  // inherits input p, h and phase
  // inherits output d
  // inherits algorithm

    annotation(derivative=density_ph_der);
  end density_ph;


  redeclare function extends density_derp_h
  "returns density derivative (dd/dp)@h=const"
  //input state
  //output ddph

protected
    HelmholtzDerivs f;
    SaturationProperties sat;
    Types.DerPressureByDensity dpdT;
    Types.DerPressureByTemperature dpTd= pressure_derT_d(state=state, f=f, sat=sat);
          DerDensityByTemperature ddTh;

  algorithm
    if (state.phase == 1) then
      f:=setHelmholtzDerivs(T=state.T, d=state.d, phase=1);
      dpdT := pressure_derd_T(state=state);
      ddTh := density_derT_h(state=state);
      ddph := 1.0/(dpdT + dpTd/ddTh);
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T, phase=2);
      ddph := state.d^2/state.T*specificHeatCapacityCv(state=state)/dpTd^2
              + state.d/state.T/dpTd;
    end if;
  end density_derp_h;


  redeclare function extends density_derh_p
  "returns density derivative (dd/dh)@p=const"
  //input state
  //output ddhp

protected
    HelmholtzDerivs f;
    SaturationProperties sat;
    Types.DerEnthalpyByDensity dhdT;
    Types.DerEnthalpyByTemperature dhTd;
          DerDensityByTemperature ddTp;

  algorithm
    if (state.phase == 1) then
      f:=setHelmholtzDerivs(
        T=state.T,
        d=state.d,
        phase=1);
      dhdT := specificEnthalpy_derd_T(state=state);
      dhTd := specificEnthalpy_derT_d(state=state);
      ddTp := density_derT_p(state=state);
      ddhp := 1.0/(dhdT + dhTd/ddTp);
    elseif (state.phase == 2) then
      sat:=setSat_T(T=state.T, phase=2);
      // dvhp = (v"-v')/(h"-h')
      // ddhp = -d^2 * dvhp
      ddhp := -state.d^2*(1/sat.liq.d-1/sat.vap.d)/(sat.liq.h-sat.vap.h);
    end if;
  end density_derh_p;


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
  // inherits input sat and output dhldp

protected
    DerEnthalpyByPressure dhpT = isothermalThrottlingCoefficient(state=sat.vap);
    DerEnthalpyByTemperature dhTp = specificHeatCapacityCp(state=sat.vap);
    DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);

  algorithm
    dhldp := dhpT + dhTp*dTp;
  end dDewEnthalpy_dPressure;
end PartialHelmholtzMedium;
