within HelmholtzFluids.Examples;
model StateRec_ps_sine
  package medium = HelmholtzFluids.R134a;
  medium.AbsolutePressure p;
  medium.SpecificEntropy s;
  medium.ThermodynamicState inletState;

protected
  constant medium.AbsolutePressure pmin=medium.fluidLimits.PMIN;
  constant medium.AbsolutePressure pcrit=medium.fluidConstants[1].criticalPressure;
  constant medium.AbsolutePressure pmax=medium.fluidLimits.PMAX;
  constant medium.SpecificEntropy smin=medium.fluidLimits.SMIN;
  constant medium.SpecificEntropy smax=medium.fluidLimits.SMAX;

Modelica.Blocks.Sources.Ramp ramp(
    height=1,
    duration=5,
    offset=0,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
Modelica.Blocks.Sources.ExpSine expSine(
    amplitude=0.99,
    offset=1,
    damping=0.5)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
Modelica.Blocks.Sources.Sine sine(
    offset=1,
    amplitude=0.999,
    phase=1.5707963267949)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

equation
  p = 0.3*pcrit + 0.3*pcrit*(ramp.y*sine.y - ramp.y);
  s = 2e3;

  inletState=medium.setState_psX(p=p, s=s, phase=0, X={1});
  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput,
    Diagram(graphics));
end StateRec_ps_sine;
