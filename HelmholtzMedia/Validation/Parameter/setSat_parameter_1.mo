within HelmholtzMedia.Validation.Parameter;

model setSat_parameter_1
  replaceable package Medium = .HelmholtzMedia.HelmholtzFluids.Propane;

  Medium.SaturationProperties sat;

  parameter Medium.Temperature T_input= 0.7*Tcrit;

protected
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;

equation
  // forward
  sat = Medium.setSat_T(T=T_input);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end setSat_parameter_1;
