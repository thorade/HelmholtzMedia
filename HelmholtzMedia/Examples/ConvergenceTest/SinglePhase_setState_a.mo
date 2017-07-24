within HelmholtzMedia.Examples.ConvergenceTest;
model SinglePhase_setState_a
  replaceable package Medium = HelmholtzFluids.Carbondioxide;
  Medium.AbsolutePressure p(start=101325);
  Medium.Temperature T(start=298.15);

  Medium.ThermodynamicState state;
  Medium.ThermodynamicState state_dT;
  Medium.ThermodynamicState state_pd;
  Medium.ThermodynamicState state_ph;
  Medium.ThermodynamicState state_ps;
  Medium.ThermodynamicState state_Ts;

Modelica.Blocks.Sources.Ramp T_sub(
    duration=4,
    startTime=0.1,
    height=Tcrit - Tmin - 3,
    offset=Tmin + 3)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

Modelica.Blocks.Sources.Ramp T_super(
    height=Tmax - Tcrit,
    duration=5,
    startTime=6,
    offset=0)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

Modelica.Blocks.Sources.Sine p_sine(
    amplitude=(pmax - pmin)/2,
    offset=(pmax - pmin)/2,
    freqHz=100,
    startTime=0)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

protected
  Medium.AbsolutePressure p_melt;
  final constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  final constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  final constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  final constant Medium.AbsolutePressure pmin=Medium.fluidConstants[1].triplePointPressure;
  final constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  final constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

equation
  T = T_sub.y + T_super.y;
  p_melt = Medium.Ancillary.meltingPressure_T(T);
  p = min(p_melt, p_sine.y);

  // set state to valid single-phase values
  state=Medium.setState_pT(p=p, T=T, phase=0);

  // call other setState functions
  state_dT=Medium.setState_dT(d=Medium.density(state), T=Medium.temperature(state), phase=0);
  state_pd=Medium.setState_pd(p=Medium.pressure(state), d=Medium.density(state), phase=0);
  state_ph=Medium.setState_ph(p=Medium.pressure(state), h=Medium.specificEnthalpy(state), phase=0);
  state_ps=Medium.setState_ps(p=Medium.pressure(state), s=Medium.specificEntropy(state), phase=0);
  state_Ts=Medium.setState_Ts(T=Medium.temperature(state), s=Medium.specificEntropy(state), phase=0);

  annotation (experiment(StopTime=12));
end SinglePhase_setState_a;
