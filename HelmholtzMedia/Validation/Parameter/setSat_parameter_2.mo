within HelmholtzMedia.Validation.Parameter;

model setSat_parameter_2
  replaceable package Medium = .HelmholtzMedia.HelmholtzFluids.Propane;

  Modelica.Units.SI.AbsolutePressure p_sat;

  parameter Medium.Temperature T_input= 0.7*Tcrit;

protected
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;

equation
  // forward
  p_sat = Medium.Internal.Ancillary.saturationPressure_T(T=T_input);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end setSat_parameter_2;
