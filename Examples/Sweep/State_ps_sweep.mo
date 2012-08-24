within HelmholtzMedia.Examples.Sweep;
model State_ps_sweep
  package Medium = HelmholtzFluids.R134a;
  Medium.AbsolutePressure p;
  Medium.SpecificEntropy s;
  Medium.ThermodynamicState inletState;

protected
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;
  constant Medium.SpecificEntropy smin=Medium.fluidLimits.SMIN;
  constant Medium.SpecificEntropy smax=Medium.fluidLimits.SMAX;

Modelica.Blocks.Sources.Ramp ramp(
    height=1,
    duration=5,
    offset=0,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
Modelica.Blocks.Sources.ExpSine expSine(
    amplitude=0.99,
    offset=1,
    damping=0.5,
    freqHz=2)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
Modelica.Blocks.Sources.Sine sine(
    offset=1,
    amplitude=0.999,
    freqHz=1,
    phase=1.5707963267949)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

equation
  p = 0.3*pcrit + 0.3*pcrit*(ramp.y*sine.y - ramp.y);
  s = 2e3;

  inletState=Medium.setState_psX(p=p, s=s, phase=0, X={1});
  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput,
    Diagram(graphics));
end State_ps_sweep;
