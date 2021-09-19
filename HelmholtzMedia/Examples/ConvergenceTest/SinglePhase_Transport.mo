within HelmholtzMedia.Examples.ConvergenceTest;
model SinglePhase_Transport
  replaceable package Medium = HelmholtzFluids.Carbondioxide;

  Medium.AbsolutePressure p;
  Medium.Temperature T;
  Medium.ThermodynamicState state;
  Medium.DynamicViscosity eta;
  Medium.ThermalConductivity lambda;

Modelica.Blocks.Sources.Ramp T_sub(
    height=Tcrit - Tmin,
    duration=4,
    offset=Tmin,
    startTime=0.1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
Modelica.Blocks.Sources.Sine p_sine(
    phase=0,
    amplitude=(pmax - pmin)/2,
    offset=(pmax - pmin)/2,
    startTime=0,
    freqHz=100)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

Modelica.Blocks.Sources.Ramp T_super(
    height=Tmax - Tcrit,
    duration=5,
    startTime=6,
    offset=0)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

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

  state=Medium.setState_pT(p=p, T=T, phase=0);
  eta=Medium.dynamicViscosity(state);
  lambda=Medium.thermalConductivity(state);

  annotation (experiment(StopTime=12));
end SinglePhase_Transport;
