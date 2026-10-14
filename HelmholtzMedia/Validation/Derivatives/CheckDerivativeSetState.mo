within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeSetState
  "time derivatives of setState_phX, setState_pTX, setState_dTX, setState_psX and of functions of a state computed by setState_phX"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
  parameter EndStates endStates=EndStates.SinglePhase "end states of the linear ramps" annotation(Evaluate=true);
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
  constant Integer n=22 "number of tested expressions";
  constant String names[n]={
    "phX.temperature","phX.density","phX.specificInternalEnergy","phX.specificEntropy",
    "density_derh_p","density_derp_h","specificHeatCapacityCp","isothermalThrottlingCoefficient",
    "pressure_derT_d","specificEnthalpy_derT_d",
    "pTX.density","pTX.specificEnthalpy","pTX.specificInternalEnergy","pTX.specificEntropy",
    "dTX.pressure","dTX.specificEnthalpy","dTX.specificInternalEnergy","dTX.specificEntropy",
    "psX.density","psX.temperature","psX.specificEnthalpy","psX.specificInternalEnergy"};
  final parameter Boolean skip_state=(endStates==EndStates.PhaseChange)
    "the derivative properties jump at the phase boundary" annotation(Evaluate=true);
  final parameter Boolean skip_inf=(endStates==EndStates.TwoPhase)
    "cp and the isothermal throttling coefficient are infinite in the two-phase region" annotation(Evaluate=true);
  final parameter Boolean skip_pT=(endStates==EndStates.TwoPhase or endStates==EndStates.PhaseChange)
    "p and T are not independent in the two-phase region" annotation(Evaluate=true);

  Real r(unit="1")=time/tEnd;
  Medium.AbsolutePressure p=s0.p + (s1.p - s0.p)*r;
  Medium.Temperature T=s0.T + (s1.T - s0.T)*r;
  Medium.SpecificEnthalpy h=s0.h + (s1.h - s0.h)*r;
  Medium.SpecificEntropy s=s0.s + (s1.s - s0.s)*r;
  Medium.Density d=s0.d + (s1.d - s0.d)*r;
  Real y[n] "expressions of setState_xxX";
  Real y_int[n] "integrated time derivatives";
  parameter Real y_start[n](each fixed=false) "absolute values at the start";
  Real err[n]={abs(y_int[i] - y[i])/(abs(y[i]) + max(1e-3*y_start[i], Modelica.Constants.small)) for i in 1:n}
    "relative deviation (offset 1e-3 of the start value)";
equation
  // no record variable for the state: its Integer phase would be a non-Real variable computed in continuous time
  y[1] = Medium.temperature(Medium.setState_phX(p=p, h=h));
  y[2] = Medium.density(Medium.setState_phX(p=p, h=h));
  y[3] = Medium.specificInternalEnergy(Medium.setState_phX(p=p, h=h));
  y[4] = Medium.specificEntropy(Medium.setState_phX(p=p, h=h));
  y[5] = if skip_state then 0 else Medium.density_derh_p(Medium.setState_phX(p=p, h=h));
  y[6] = if skip_state then 0 else Medium.density_derp_h(Medium.setState_phX(p=p, h=h));
  y[7] = if skip_state or skip_inf then 0 else Medium.specificHeatCapacityCp(Medium.setState_phX(p=p, h=h));
  y[8] = if skip_state or skip_inf then 0 else Medium.isothermalThrottlingCoefficient(Medium.setState_phX(p=p, h=h));
  y[9] = if skip_state then 0 else Medium.pressure_derT_d(Medium.setState_phX(p=p, h=h));
  y[10] = if skip_state then 0 else Medium.specificEnthalpy_derT_d(Medium.setState_phX(p=p, h=h));
  y[11] = if skip_pT then 0 else Medium.density(Medium.setState_pTX(p=p, T=T));
  y[12] = if skip_pT then 0 else Medium.specificEnthalpy(Medium.setState_pTX(p=p, T=T));
  y[13] = if skip_pT then 0 else Medium.specificInternalEnergy(Medium.setState_pTX(p=p, T=T));
  y[14] = if skip_pT then 0 else Medium.specificEntropy(Medium.setState_pTX(p=p, T=T));
  y[15] = Medium.pressure(Medium.setState_dTX(d=d, T=T));
  y[16] = Medium.specificEnthalpy(Medium.setState_dTX(d=d, T=T));
  y[17] = Medium.specificInternalEnergy(Medium.setState_dTX(d=d, T=T));
  y[18] = Medium.specificEntropy(Medium.setState_dTX(d=d, T=T));
  y[19] = Medium.density(Medium.setState_psX(p=p, s=s));
  y[20] = Medium.temperature(Medium.setState_psX(p=p, s=s));
  y[21] = Medium.specificEnthalpy(Medium.setState_psX(p=p, s=s));
  y[22] = Medium.specificInternalEnergy(Medium.setState_psX(p=p, s=s));
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
Tests the derivative annotations of setState_phX, setState_pTX, setState_dTX and setState_psX (Internal.setState_*_der)
and, through them, the functions of a state computed by setState_phX, as in Modelica.Fluid components. All inputs vary linearly between the end states selected by <code>endStates</code>
(see CheckDerivativeAnnotations).
</p>
</html>"));
end CheckDerivativeSetState;
