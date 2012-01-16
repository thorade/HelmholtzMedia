within HelmholtzFluids.Examples;
model StateRec_ph_sine
  package medium = HelmholtzFluids.Butane;
  medium.AbsolutePressure p;
  medium.SpecificEnthalpy h;
  medium.ThermodynamicState inletState;

  Modelica.Blocks.Sources.Sine sine_p(
    freqHz=1.1,
    startTime=5,
    amplitude=1.9e6,
    offset=2e6,
    phase=0)
    annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
  Modelica.Blocks.Sources.ExpSine expSine_h(
    damping=0.3,
    amplitude=500e3,
    offset=-200e3,
    startTime=0,
    freqHz=0.8,
    phase=1.5707963267949)
    annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));

equation
  p = sine_p.y;
  h = expSine_h.y;
  inletState=medium.setState_phX(p=p, h=h, phase=0, X={1});
  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end StateRec_ph_sine;
