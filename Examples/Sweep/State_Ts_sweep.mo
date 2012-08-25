within HelmholtzMedia.Examples.Sweep;
model State_Ts_sweep
  package Medium = HelmholtzFluids.R134a;
  Medium.Temperature T;
  Medium.SpecificEntropy s;
  Medium.ThermodynamicState inletState;

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  constant Medium.SpecificEntropy smin=Medium.fluidLimits.SMIN;
  constant Medium.SpecificEntropy smax=Medium.fluidLimits.SMAX;

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
  T = 0.7*Tcrit + (expSine.y- 0.1*ramp.y)*(Tcrit-Tmin);
  s = 2e3 + 1e2*(ramp.y*sine.y - ramp.y);

  inletState=Medium.setState_Ts(T=T, s=s, phase=0);
  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput,
    Diagram(graphics));
end State_Ts_sweep;
