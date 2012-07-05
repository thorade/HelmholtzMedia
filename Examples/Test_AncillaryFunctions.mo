within HelmholtzMedia.Examples;
model Test_AncillaryFunctions
  package Medium = HelmholtzFluids.Ethanol;
  Medium.Temperature Tsat;
  Medium.AbsolutePressure psat;
  Medium.Density dliq;
  Medium.Density dvap;

protected
  constant Medium.Temperature T_trip=Medium.fluidConstants[1].triplePointTemperature;
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;

  Modelica.Blocks.Sources.Ramp Temp_ramp(
    height=T_crit-T_trip-(1e-12),
    duration=5,
    offset=T_trip,
    startTime=2)
    annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

algorithm
    Tsat := Temp_ramp.y;
    psat := Medium.saturationPressure(T=Tsat);
    dliq := Medium.bubbleDensity_T_ANC(T=Tsat);
    dvap := Medium.dewDensity_T_ANC(T=Tsat);

  annotation (experiment(StopTime=10), __Dymola_experimentSetupOutput);
end Test_AncillaryFunctions;
