within HelmholtzFluids.Examples;
model StateRec_dT_sine
  package medium = HelmholtzFluids.Butane;
  medium.Density d;
  medium.Temperature T;
  medium.ThermodynamicState inletState;
  medium.AbsolutePressure p_in;
  medium.SpecificEnthalpy h_in;

  Modelica.Blocks.Sources.ExpSine expSine_d(
    offset=300,
    damping=0.5,
    startTime=3,
    amplitude=299,
    freqHz=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine sine_T(
    startTime=1,
    freqHz=4,
    offset=298.15,
    amplitude=150)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

equation
  d = expSine_d.y;
  T = sine_T.y;
  inletState=medium.setState_dTX(d=d, T=T, phase=0, X={1});
  p_in=inletState.p;
  h_in=inletState.h;

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end StateRec_dT_sine;
