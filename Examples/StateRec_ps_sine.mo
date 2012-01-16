within HelmholtzFluids.Examples;
model StateRec_ps_sine
  package medium = HelmholtzFluids.Butane;
  medium.AbsolutePressure p;
  medium.SpecificEntropy s;
  medium.ThermodynamicState inletState;

  Modelica.Blocks.Sources.Sine sine_p(
    freqHz=1.1,
    startTime=5,
    phase=0,
    amplitude=0.9e6,
    offset=1e6)
    annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
  Modelica.Blocks.Sources.ExpSine expSine_s(
    damping=0.3,
    startTime=0,
    freqHz=0.8,
    offset=0,
    amplitude=1e3,
    phase=1.5707963267949)
    annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));

equation
  p = 2e6; //sine_p.y;
  s = expSine_s.y;
  inletState=medium.setState_psX(p=p, s=s, phase=0, X={1});
  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput,
    Diagram(graphics));
end StateRec_ps_sine;
