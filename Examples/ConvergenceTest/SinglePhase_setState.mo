within HelmholtzMedia.Examples.ConvergenceTest;
model SinglePhase_setState
  package Medium = HelmholtzFluids.Pentane;
  Medium.AbsolutePressure p;
  Medium.Temperature T;
  Medium.ThermodynamicState state;
  Medium.ThermodynamicState state_dT;
  Medium.ThermodynamicState state_pd;
  Medium.ThermodynamicState state_ph;
  Medium.ThermodynamicState state_ps;
  Medium.ThermodynamicState state_Ts;

Modelica.Blocks.Sources.Ramp T_sub(
    duration=4,
    startTime=0.1,
    height=Tcrit - Tmin,
    offset=Tmin)
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
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

equation
  T = T_sub.y + T_super.y;
  p = p_sine.y;

  // set state to valid single-phase values
  state=Medium.setState_pT(p=p, T=T, phase=0);

  // call other setState functions
  state_dT=Medium.setState_dT(d=state.d, T=state.T, phase=0);
  state_pd=Medium.setState_pd(p=state.p, d=state.d, phase=0);
  state_ph=Medium.setState_ph(p=state.p, h=state.h, phase=0);
  state_ps=Medium.setState_ps(p=state.p, s=state.s, phase=0);
  state_Ts=Medium.setState_Ts(T=state.T, s=state.s, phase=0);

  annotation (experiment(StopTime=12, NumberOfIntervals=10000));
end SinglePhase_setState;
