within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotations
  "time derivatives of the property functions of two variables (derivative annotations) against the functions themselves"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
  parameter EndStates endStates=EndStates.SinglePhase "end states of the linear ramps" annotation(Evaluate=true);
  final parameter Boolean skip_pT=(endStates==EndStates.TwoPhase or endStates==EndStates.PhaseChange)
    "the functions of p and T are skipped if the path enters the two-phase region (p and T not independent there)" annotation(Evaluate=true);
  parameter Real tol=if endStates==EndStates.PhaseChange then 1e-5 else 1e-6
    "relative tolerance (larger for PhaseChange: the derivatives jump at the phase boundaries, without events)";
  parameter Modelica.Units.SI.Time tEnd=1 "duration of the ramps";
  parameter Real dT_sat=1e-4 "relative distance of the end state to the saturation temperature (LiquidNearSaturation, VapourNearSaturation)";
  parameter Medium.SaturationProperties sat0=Medium.setSat_p(p=0.5*Medium.fluidConstants[1].criticalPressure);
  parameter Medium.SaturationProperties sat1=Medium.setSat_p(p=0.6*Medium.fluidConstants[1].criticalPressure);
  parameter Medium.ThermodynamicState s0=
    if endStates==EndStates.TwoPhase then Medium.setState_ph(p=sat0.psat, h=sat0.liq.h + 0.3*(sat0.vap.h - sat0.liq.h))
    elseif endStates==EndStates.PhaseChange then Medium.setState_pT(p=sat0.psat, T=0.95*sat0.Tsat)
    elseif endStates==EndStates.LiquidNearSaturation then Medium.setState_pT(p=sat0.psat, T=0.9*sat0.Tsat)
    elseif endStates==EndStates.VapourNearSaturation then Medium.setState_pT(p=sat0.psat, T=(1 + dT_sat)*sat0.Tsat)
    else Medium.setState_pT(p=1.1*Medium.fluidConstants[1].criticalPressure, T=1.05*Medium.fluidConstants[1].criticalTemperature);
  parameter Medium.ThermodynamicState s1=
    if endStates==EndStates.TwoPhase then Medium.setState_ph(p=sat1.psat, h=sat1.liq.h + 0.7*(sat1.vap.h - sat1.liq.h))
    elseif endStates==EndStates.PhaseChange then Medium.setState_pT(p=sat0.psat, T=1.05*sat0.Tsat)
    elseif endStates==EndStates.LiquidNearSaturation then Medium.setState_pT(p=sat0.psat, T=(1 - dT_sat)*sat0.Tsat)
    elseif endStates==EndStates.VapourNearSaturation then Medium.setState_pT(p=sat0.psat, T=1.1*sat0.Tsat)
    else Medium.setState_pT(p=1.5*Medium.fluidConstants[1].criticalPressure, T=1.3*Medium.fluidConstants[1].criticalTemperature);
  constant String names[22]={"density_ph","temperature_ph","specificEntropy_ph","density_ps","temperature_ps","specificEnthalpy_ps",
    "pressure_dT","specificEnthalpy_dT","specificEntropy_dT","pressure_Ts","density_Ts","temperature_pd","specificEnthalpy_pd",
    "density_pT","specificEnthalpy_pT","specificEntropy_pT",
    "density_derh_p","density_derp_h","specificHeatCapacityCp","isothermalThrottlingCoefficient",
    "dynamicViscosity","thermalConductivity"};
  constant Integer n=22 "number of tested functions";
  final parameter Boolean skip_state=(endStates==EndStates.PhaseChange)
    "the functions of the state are skipped if the path crosses the phase boundary (fixed phase of the state)" annotation(Evaluate=true);
  final parameter Medium.FixedPhase phase_state=if endStates==EndStates.TwoPhase then 2 else 1
    "phase of the state for the functions of the state" annotation(Evaluate=true);

  // each pair of inputs varies linearly between the two end states
  Real r(unit="1")=time/tEnd;
  Medium.AbsolutePressure p=s0.p + (s1.p - s0.p)*r;
  Medium.Temperature T=s0.T + (s1.T - s0.T)*r;
  Medium.SpecificEnthalpy h=s0.h + (s1.h - s0.h)*r;
  Medium.SpecificEntropy s=s0.s + (s1.s - s0.s)*r;
  Medium.Density d=s0.d + (s1.d - s0.d)*r;
  // state built from variables (as in BaseProperties), so that der() of the functions of the state is defined
  Medium.Density d_state=Medium.density_ph(p=p, h=h);
  Medium.Temperature T_state=Medium.temperature_ph(p=p, h=h);
  Medium.ThermodynamicState state=Medium.ThermodynamicState(
    phase=phase_state, d=d_state, T=T_state, p=p, h=h, u=h - p/d_state, s=Medium.specificEntropy_ph(p=p, h=h));
  Real y[n] "property functions";
  Real y_int[n] "integrated time derivatives";
  parameter Real y_start[n](each fixed=false) "absolute values at the start, scale of the functions of the state";
  Real err[n]={abs(y_int[i] - y[i])/(abs(y[i]) + (if i <= 16 then 1e-3 else max(1e-3*y_start[i], Modelica.Constants.small))) for i in 1:n}
    "relative deviation (offset 1e-3 for the functions of two variables, 1e-3 of the start value for the functions of the state)";
equation
  y[1] = Medium.density_ph(p=p, h=h);
  y[2] = Medium.temperature_ph(p=p, h=h);
  y[3] = Medium.specificEntropy_ph(p=p, h=h);
  y[4] = Medium.density_ps(p=p, s=s);
  y[5] = Medium.temperature_ps(p=p, s=s);
  y[6] = Medium.specificEnthalpy_ps(p=p, s=s);
  y[7] = Medium.pressure_dT(d=d, T=T);
  y[8] = Medium.specificEnthalpy_dT(d=d, T=T);
  y[9] = Medium.specificEntropy_dT(d=d, T=T);
  y[10] = Medium.pressure_Ts(T=T, s=s);
  y[11] = Medium.density_Ts(T=T, s=s);
  y[12] = Medium.temperature_pd(p=p, d=d);
  y[13] = Medium.specificEnthalpy_pd(p=p, d=d);
  y[14] = if skip_pT then 0 else Medium.density_pT(p=p, T=T);
  y[15] = if skip_pT then 0 else Medium.specificEnthalpy_pT(p=p, T=T);
  y[16] = if skip_pT then 0 else Medium.specificEntropy_pT(p=p, T=T);
  y[17] = if skip_state then 0 else Medium.density_derh_p(state);
  y[18] = if skip_state then 0 else Medium.density_derp_h(state);
  y[19] = if skip_state or phase_state==2 then 0 else Medium.specificHeatCapacityCp(state) "infinite in the two-phase region";
  y[20] = if skip_state or phase_state==2 then 0 else Medium.isothermalThrottlingCoefficient(state) "infinite in the two-phase region";
  y[21] = if skip_state or phase_state==2 then 0 else Medium.dynamicViscosity(state) "not defined in the two-phase region";
  y[22] = if skip_state or phase_state==2 then 0 else Medium.thermalConductivity(state) "not defined in the two-phase region";
  der(y_int) = der(y);
initial equation
  y_int = y;
  y_start = abs(y);
algorithm
  when terminal() then
    for i in 1:n loop
      Modelica.Utilities.Streams.print(getInstanceName() + "." + names[i] + ": relative deviation " + String(err[i], significantDigits=3));
    end for;
    assert(max(err) < tol, getInstanceName() + ": a time derivative deviates, maximum relative deviation " + String(max(err)));
  end when;
  annotation (experiment(StopTime=1, Tolerance=1e-10), Documentation(info="<html>
<p>
For each property function of two variables (density_ph, pressure_Ts, ...) the time derivative computed by the tool
from the derivative annotation (functions *_state and *_der in PartialHelmholtzMedium.Internal) is integrated and
compared to the function value. The deviation is relative to the value of each property (plus 1e-3).
</p>
<p>
All inputs vary linearly between the end states s0 and s1, selected by <code>endStates</code>:
</p>
<ul>
<li>SinglePhase: supercritical, 1.1 p_crit, 1.05 T_crit to 1.5 p_crit, 1.3 T_crit.</li>
<li>TwoPhase: two-phase states at 0.5 p_crit (vapour quality 0.3) and 0.6 p_crit (vapour quality 0.7).</li>
<li>PhaseChange: at 0.5 p_crit from subcooled liquid (0.95 T_sat) through the two-phase region to superheated
vapour (1.05 T_sat); the derivatives jump at the phase boundaries.</li>
<li>LiquidNearSaturation: at 0.5 p_crit from 0.9 T_sat to (1 - dT_sat) T_sat.</li>
<li>VapourNearSaturation: at 0.5 p_crit from (1 + dT_sat) T_sat to 1.1 T_sat.</li>
</ul>
<p>
The functions of the state (density_derh_p, density_derp_h, specificHeatCapacityCp, isothermalThrottlingCoefficient,
dynamicViscosity, thermalConductivity) are evaluated for a state built from variables, as in BaseProperties, with a fixed phase; they are skipped for
PhaseChange, cp, the isothermal throttling coefficient and the transport properties also for TwoPhase (infinite or not defined there).
Their deviation is relative to 1e-3 of their start value, since their values are small in SI units.
</p>
<p>
The functions of p and T are skipped for TwoPhase and PhaseChange, since p and T are not independent in the
two-phase region (density_pT etc. jump at the saturation line). The other paths may touch the two-phase region,
the functions of the other variable pairs are continuous there.
</p>
</html>"));
end CheckDerivativeAnnotations;
