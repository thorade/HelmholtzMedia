within HelmholtzFluids.Examples;
model StateRec_pT_sine
  package medium = HelmholtzFluids.Butane;
  medium.AbsolutePressure p;
  medium.Temperature T;
  medium.ThermodynamicState inletState;
  medium.Density d_in;
  medium.SpecificEnthalpy h_in;
  medium.DynamicViscosity eta;
  medium.ThermalConductivity lambda;
  medium.PrandtlNumber Pr;

  Modelica.Blocks.Sources.ExpSine expSine_p(
    offset=2e6,
    damping=0.5,
    startTime=3,
    amplitude=1.99e6,
    freqHz=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine sine_T(
    startTime=1,
    freqHz=4,
    amplitude=150,
    offset=298.15)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

equation
  p = expSine_p.y;
  T = sine_T.y;
  inletState=medium.setState_pTX(p=p, T=T, phase=0, X={1});
  d_in=inletState.d;
  h_in=inletState.h;
  eta=medium.dynamicViscosity(inletState);
  lambda=medium.thermalConductivity(inletState);
  Pr=medium.prandtlNumber(inletState);

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end StateRec_pT_sine;
