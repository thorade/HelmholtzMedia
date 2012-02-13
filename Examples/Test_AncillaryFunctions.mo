within HelmholtzFluids.Examples;
model Test_AncillaryFunctions
  package medium = HelmholtzFluids.Butane;
  medium.Temperature Tsat;
  medium.AbsolutePressure psat;
  medium.Density dliq;
  medium.Density dvap;
  constant medium.Temperature T_trip=medium.fluidConstants[1].triplePointTemperature;
  constant medium.Temperature T_crit=medium.fluidConstants[1].criticalTemperature;

  Modelica.Blocks.Sources.Ramp Temp_ramp(
    height=T_crit-T_trip,
    duration=8,
    offset=T_trip,
    startTime=1)
    annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

algorithm
    Tsat := Temp_ramp.y;
    psat := medium.saturationPressure(T=Tsat);
    dliq := medium.bubbleDensity_T_ANC(T=Tsat);
    dvap := medium.dewDensity_T_ANC(T=Tsat);

  annotation (experiment(StopTime=10), __Dymola_experimentSetupOutput);
end Test_AncillaryFunctions;
