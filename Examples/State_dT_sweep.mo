within HelmholtzMedia.Examples;
model State_dT_sweep
  package Medium = HelmholtzFluids.R134a;
  Medium.Density d;
  Medium.Temperature T;
  Medium.ThermodynamicState inletState;
  Medium.AbsolutePressure p_in;
  Medium.SpecificEnthalpy h_in;

protected
  constant Medium.Density dmin=Medium.fluidLimits.DMIN;
  constant Medium.Density dcrit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Density dmax=Medium.fluidLimits.DMAX;
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;

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
  d = dcrit - ramp.y*dcrit + ramp.y*sine.y*dcrit;
  T = 0.7*Tcrit + expSine.y*(Tcrit-Tmin);

  inletState=Medium.setState_dTX(d=d, T=T, phase=0, X={1});
  p_in=inletState.p;
  h_in=inletState.h;

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end State_dT_sweep;
