within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_Ts "Return thermodynamic state as function of (T, s)"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Entropy";
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
    assert(T >= T_trip, "setState_Ts_error: pressure is lower than triple point pressure");
    assert(T <= T_crit, "setState_TsX_error: pressure is higher than critical pressure");
    sat := setSat_T(T=T);
    assert(s >= sat.liq.s, "setState_TsX_error: entropy is lower than saturated liquid entropy: this is single phase liquid");
    assert(s <= sat.vap.s, "setState_TsX_error: entropy is higher than saturated vapor entropy: this is single phase vapor");
  else
    if ((T <= T_crit) and (T >= T_trip)) then
      // two-phase possible, do simple check first
      sat.Tsat := T;
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
        // two-phase state or close to it, get saturation properties from EoS
        sat := setSat_T(T=sat.Tsat);
      end if;

      if (s < sat.liq.s) then
        state.phase := 1; // single phase liquid
        dmin := sat.liq.d;
      elseif (s > sat.vap.s) then
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
  state.s := s;
  if (state.phase == 2) then
    // force two-phase, SaturationProperties are already known
    state.p := sat.psat;
    x := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
    state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
  else
    // force single-phase
    state.d := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function setState_Ts_RES(
            T=T,
            s=s,
            phase=1),
          u_min=0.98*dmin,
          u_max=1.02*dmax,
          tolerance=tolerance);

    tau := T_crit/state.T;
    delta := state.d/d_crit;

    f.it  := EoS.f_it(tau=tau, delta=delta);
    f.rt  := EoS.f_rt(tau=tau, delta=delta);
    f.rd  := EoS.f_rd(tau=tau, delta=delta);
    state.p := state.d*T*R*(1+delta*f.rd);
    state.h := state.T*R*(tau*(f.it + f.rt) + (1+delta*f.rd));
    state.u := state.T*R*(tau*(f.it+f.rt));
  end if;

end setState_Ts;
