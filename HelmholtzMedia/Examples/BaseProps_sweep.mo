within HelmholtzMedia.Examples;
model BaseProps_sweep
  "calculate BaseProperties from any two given input properties"

  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  // replaceable package Medium = HelmholtzFluids.R134a(dT_explicit=true);
  // replaceable package Medium = HelmholtzFluids.Butane(independentVariables=IndependentVariables.ph);
  // replaceable package Medium = HelmholtzFluids.R134a(independentVariables=IndependentVariables.ph);
  Medium.BaseProperties medium;

  Modelica.Blocks.Sources.ExpSine d_generator(
    amplitude=50,
    offset=300,
    damping=0.5,
    f=1.1,
    startTime=2) annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine T_generator(
    offset=298.15,
    f=1.3,
    startTime=0,
    amplitude=50) annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Sine p_generator(
    offset=2e6,
    startTime=0,
    phase=0,
    f=1.3,
    amplitude=1.999e6) annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
  Modelica.Blocks.Sources.ExpSine h_generator(
    damping=0.3,
    f=2,
    amplitude=100e3,
    startTime=2,
    phase=1.5707963267949,
    offset=500e3) annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));

  Modelica.Blocks.Sources.ExpSine s_generator(
    damping=0.3,
    startTime=2,
    phase=0,
    f=1.1,
    amplitude=2e3,
    offset=0) annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));

equation
  // medium.d = d_generator.y;
  // medium.T = T_generator.y;
  // medium.s = s_generator.y;
  medium.p = p_generator.y;
  medium.h = h_generator.y;

  annotation (experiment(
      StopTime=10,
      Tolerance=1e-05));
end BaseProps_sweep;
