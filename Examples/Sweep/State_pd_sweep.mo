within HelmholtzMedia.Examples.Sweep;
model State_pd_sweep
  package Medium = HelmholtzFluids.Butane;
  Medium.AbsolutePressure p;
  Medium.Density d;
  Medium.ThermodynamicState state;

protected
  constant Medium.Density dmin=Medium.fluidLimits.DMIN;
  constant Medium.Density dcrit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Density dmax=Medium.fluidLimits.DMAX;
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

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
  p = pcrit + pcrit*(ramp.y*sine.y - ramp.y);
  d = dcrit - ramp.y*dcrit + ramp.y*sine.y*dcrit;
  state=Medium.setState_pd(p=p, d=d, phase=0);

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end State_pd_sweep;
