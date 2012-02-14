within HelmholtzFluids.Examples;
model BaseProps_sweep
  "calculate BaseProperties from any two given input properties"
  package medium = HelmholtzFluids.Butane;
  medium.BaseProperties props;

  Modelica.Blocks.Sources.ExpSine d_generator(
    amplitude=50,
    offset=300,
    damping=0.5,
    freqHz=1.1,
    startTime=2)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine T_generator(
    amplitude=125,
    offset=298.15,
    freqHz=1.3,
    startTime=0)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Sine p_generator(
    offset=2e6,
    startTime=0,
    phase=0,
    freqHz=1.3,
    amplitude=1.999e6)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
  Modelica.Blocks.Sources.ExpSine h_generator(
    damping=0.3,
    offset=-200e3,
    freqHz=2,
    amplitude=100e3,
    phase=1.5707963267949,
    startTime=2)
    annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));

  Modelica.Blocks.Sources.ExpSine s_generator(
    damping=0.3,
    startTime=2,
    phase=0,
    freqHz=1.1,
    amplitude=2e3,
    offset=0)
    annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));
equation
  // props.d = d_generator.y;
  // props.T = T_generator.y;
  // props.s = s_generator.y;
  props.p = p_generator.y;
  props.h = h_generator.y;

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput,
    Diagram(graphics));
end BaseProps_sweep;
