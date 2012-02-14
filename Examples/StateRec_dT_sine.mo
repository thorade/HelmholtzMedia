within HelmholtzFluids.Examples;
model StateRec_dT_sine
  package medium = HelmholtzFluids.R134a;
  medium.Density d;
  medium.Temperature T;
  medium.ThermodynamicState inletState;
  medium.AbsolutePressure p_in;
  medium.SpecificEnthalpy h_in;

protected
  constant medium.Density dmin=medium.fluidLimits.DMIN;
  constant medium.Density dcrit=medium.fluidConstants[1].molarMass/medium.fluidConstants[1].criticalMolarVolume;
  constant medium.Density dmax=medium.fluidLimits.DMAX;
  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature Tmax=medium.fluidLimits.TMAX;

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
  d = dcrit - ramp.y*dcrit + ramp.y*sine.y*dcrit;
  T = 0.7*Tcrit + expSine.y*(Tcrit-Tmin);

  inletState=medium.setState_dTX(d=d, T=T, phase=0, X={1});
  p_in=inletState.p;
  h_in=inletState.h;

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end StateRec_dT_sine;
