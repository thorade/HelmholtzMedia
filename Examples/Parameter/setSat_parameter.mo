within HelmholtzMedia.Examples.Parameter;
model setSat_parameter
  package medium = HelmholtzFluids.Butane;

  medium.SaturationProperties sat;
  medium.Temperature T_test;

  parameter Real T_input = 423.15;
  parameter medium.AbsolutePressure p_input = 101325;

equation
sat = medium.setSat_T(T=T_input);
// sat = medium.setSat_p(p=p_input);
T_test = medium.Ancillary.saturationTemperature_d(sat.liq.d);

  annotation (experiment(Tolerance=1e-012), __Dymola_experimentSetupOutput);
end setSat_parameter;
