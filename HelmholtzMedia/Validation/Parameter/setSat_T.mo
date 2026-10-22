within HelmholtzMedia.Validation.Parameter;
model setSat_T
  replaceable package Medium = .HelmholtzMedia.HelmholtzFluids.Propane;

  Medium.SaturationProperties sat;

  parameter Medium.Temperature T_input= 273.15;

equation
  // forward
  sat = Medium.setSat_T(T=T_input);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end setSat_T;
