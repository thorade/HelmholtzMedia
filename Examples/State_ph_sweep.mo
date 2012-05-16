within HelmholtzMedia.Examples;
model State_ph_sweep
  package Medium = HelmholtzFluids.R134a;
  Medium.AbsolutePressure p;
  Medium.SpecificEnthalpy h;
  Medium.ThermodynamicState inletState;

protected
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;
  constant Medium.SpecificEnthalpy hmin=Medium.fluidLimits.HMIN;
  constant Medium.SpecificEnthalpy hmax=Medium.fluidLimits.HMAX;

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
  p = 0.5*pcrit + 0.5*pcrit*(ramp.y*sine.y - ramp.y);
  h = 200e3;// (hmax-hmin)/2 + 0.002*expSine.y*(hmax-hmin);

  inletState=Medium.setState_phX(p=p, h=h, phase=0, X={1});
  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end State_ph_sweep;
