within HelmholtzMedia.Examples.Sweep;
model BaseProps_sweep
  "calculate BaseProperties from any two given input properties"
  //package Medium = HelmholtzFluids.Butane;
  package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  //package Medium = HelmholtzFluids.R134a(dT_explicit=true);
  //package Medium = HelmholtzFluids.Butane(independentVariables=IndependentVariables.ph);
  //package Medium = HelmholtzFluids.R134a(independentVariables=IndependentVariables.ph);
  Medium.BaseProperties props;

  Modelica.Blocks.Sources.ExpSine d_generator(
    amplitude=50,
    offset=300,
    damping=0.5,
    freqHz=1.1,
    startTime=2)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine T_generator(
    offset=298.15,
    freqHz=1.3,
    startTime=0,
    amplitude=50)
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
    freqHz=2,
    amplitude=100e3,
    startTime=2,
    phase=1.5707963267949,
    offset=500e3)
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

end BaseProps_sweep;
