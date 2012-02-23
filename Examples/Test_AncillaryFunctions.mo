within HelmholtzFluids.Examples;
model Test_AncillaryFunctions
  package medium = HelmholtzFluids.Isobutane;
  medium.Temperature Tsat;
  medium.AbsolutePressure psat;
  medium.Density dliq;
  medium.Density dvap;

protected
  constant medium.Temperature T_trip=medium.fluidConstants[1].triplePointTemperature;
  constant medium.Temperature T_crit=medium.fluidConstants[1].criticalTemperature;

  Modelica.Blocks.Sources.Ramp Temp_ramp(
    height=T_crit-T_trip-(1e-12),
    duration=5,
    offset=T_trip,
    startTime=2)
    annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

algorithm
    Tsat := Temp_ramp.y;
    psat := medium.saturationPressure(T=Tsat);
    dliq := medium.bubbleDensity_T_ANC(T=Tsat);
    dvap := medium.dewDensity_T_ANC(T=Tsat);

end Test_AncillaryFunctions;
