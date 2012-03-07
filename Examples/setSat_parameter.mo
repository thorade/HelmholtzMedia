within HelmholtzMedia.Examples;
model setSat_parameter
  package medium = HelmholtzFluids.R134a;

  medium.SaturationProperties satProps;

  parameter medium.Temperature T_input = 273.15;
  parameter medium.AbsolutePressure p_input = 101325;

equation
  satProps = medium.setSat_T(T=T_input);
//satProps = medium.setSat_p(p=p_input);

  annotation (experiment(Tolerance=1e-012), __Dymola_experimentSetupOutput);
end setSat_parameter;
