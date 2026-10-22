within HelmholtzMedia.Validation.Parameter;
model setSat_p_internal
  replaceable package Medium = .HelmholtzMedia.HelmholtzFluids.Propane;

  Modelica.Units.SI.AbsolutePressure p_sat;

  parameter Medium.Temperature T_input= 273.15;

equation
  // forward
  p_sat = Medium.Internal.Ancillary.saturationPressure_T(T=T_input);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end setSat_p_internal;
