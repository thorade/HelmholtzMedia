within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotations
  "time derivatives of the property functions of two variables (derivative annotations) against the functions themselves"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
  parameter EndStates endStates=EndStates.SinglePhase "end states of the linear ramps" annotation(Evaluate=true);
  final parameter Boolean skip_pT=(endStates==EndStates.TwoPhase or endStates==EndStates.PhaseChange)
    "the functions of p and T are skipped if the path enters the two-phase region (p and T not independent there)" annotation(Evaluate=true);
  parameter Real tol=1e-6 "relative tolerance";
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
  constant String names[16]={"density_ph","temperature_ph","specificEntropy_ph","density_ps","temperature_ps","specificEnthalpy_ps",
    "pressure_dT","specificEnthalpy_dT","specificEntropy_dT","pressure_Ts","density_Ts","temperature_pd","specificEnthalpy_pd",
    "density_pT","specificEnthalpy_pT","specificEntropy_pT"};

  // each pair of inputs varies linearly between the two end states
  Real r(unit="1")=time/tEnd;
  Medium.AbsolutePressure p=s0.p + (s1.p - s0.p)*r;
  Medium.Temperature T=s0.T + (s1.T - s0.T)*r;
  Medium.SpecificEnthalpy h=s0.h + (s1.h - s0.h)*r;
  Medium.SpecificEntropy s=s0.s + (s1.s - s0.s)*r;
  Medium.Density d=s0.d + (s1.d - s0.d)*r;
  Real y[16] "property functions";
  Real y_int[16] "integrated time derivatives";
  Real err[16]={abs(y_int[i] - y[i])/(abs(y[i]) + 1e-3) for i in 1:16} "relative deviation";
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
  der(y_int) = der(y);
initial equation
  y_int = y;
algorithm
  when terminal() then
    for i in 1:16 loop
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
The functions of p and T are skipped for TwoPhase and PhaseChange, since p and T are not independent in the
two-phase region (density_pT etc. jump at the saturation line). The other paths may touch the two-phase region,
the functions of the other variable pairs are continuous there.
</p>
</html>"));
end CheckDerivativeAnnotations;
