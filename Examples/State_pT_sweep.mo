within HelmholtzMedia.Examples;
model State_pT_sweep
  package medium = HelmholtzFluids.R134a;
  medium.AbsolutePressure p;
  medium.Temperature T;
  medium.ThermodynamicState inletState;
  medium.Density d_in;
  medium.SpecificEnthalpy h_in;
  medium.DynamicViscosity eta;
  medium.ThermalConductivity lambda;
  medium.PrandtlNumber Pr;

protected
  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature Tmax=medium.fluidLimits.TMAX;
  constant medium.AbsolutePressure pmin=medium.fluidLimits.PMIN;
  constant medium.AbsolutePressure pcrit=medium.fluidConstants[1].criticalPressure;
  constant medium.AbsolutePressure pmax=medium.fluidLimits.PMAX;

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
  p = pcrit + pcrit*(ramp.y*sine.y - ramp.y);
  T = 0.7*Tcrit + expSine.y*(Tcrit-Tmin);

  inletState=medium.setState_pTX(p=p, T=T, phase=0, X={1});
  d_in=inletState.d;
  h_in=inletState.h;
  eta=medium.dynamicViscosity(inletState);
  lambda=medium.thermalConductivity(inletState);
  Pr=medium.prandtlNumber(inletState);

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end State_pT_sweep;
