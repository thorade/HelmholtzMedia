within HelmholtzFluids.Examples;
model setSat_sweep
  package medium = HelmholtzFluids.R134a;
  medium.SaturationProperties satProps;

protected
  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature Tmax=medium.fluidLimits.TMAX;
  constant medium.AbsolutePressure pmin=medium.fluidLimits.PMIN;
  constant medium.AbsolutePressure pcrit=medium.fluidConstants[1].criticalPressure;
  constant medium.AbsolutePressure pmax=medium.fluidLimits.PMAX;

  Modelica.Blocks.Sources.Ramp Tramp(
    height=Tcrit - Tmin-(1e-12),
    offset=Tmin,
    duration=8,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Ramp pramp(
    height=pcrit - pmin,
    offset=pmin,
    duration=8,
    startTime=1) annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

equation
  satProps = medium.setSat_T(T=Tramp.y);
//satProps = medium.setSat_p(p=pramp.y);

  annotation (experiment(
      StopTime=10,
      NumberOfIntervals=10000,
      Tolerance=1e-008),               __Dymola_experimentSetupOutput);
end setSat_sweep;
