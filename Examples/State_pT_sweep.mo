within HelmholtzMedia.Examples;
model State_pT_sweep
  package Medium = HelmholtzFluids.R134a;
  Medium.AbsolutePressure p;
  Medium.Temperature T;
  Medium.ThermodynamicState inletState;
  Medium.Density d_in;
  Medium.SpecificEnthalpy h_in;
  Medium.DynamicViscosity eta;
  Medium.ThermalConductivity lambda;
  Medium.PrandtlNumber Pr;

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

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

  inletState=Medium.setState_pTX(p=p, T=T, phase=0, X={1});
  d_in=inletState.d;
  h_in=inletState.h;
  eta=Medium.dynamicViscosity(inletState);
  lambda=Medium.thermalConductivity(inletState);
  Pr=Medium.prandtlNumber(inletState);

  annotation (experiment(StopTime=10, NumberOfIntervals=1000),
      __Dymola_experimentSetupOutput);
end State_pT_sweep;
