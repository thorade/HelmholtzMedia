within HelmholtzMedia.Examples.Parameter;
model setSat_parameter
  package medium = HelmholtzFluids.Butane;

  medium.SaturationProperties sat;
  medium.SaturationProperties sat_p;
  medium.SaturationProperties sat_dl;
  medium.SaturationProperties sat_dv;

  parameter medium.Temperature T_input = 135;

equation
  // forward
  sat = medium.setSat_T(T=T_input);

  // backward
  sat_p  = medium.setSat_p(p=sat.psat);
  sat_dl = medium.setSat_d(d=sat.liq.d);
  sat_dv = medium.setSat_d(d=sat.vap.d);

  annotation (experiment(Tolerance=1e-012), __Dymola_experimentSetupOutput);
end setSat_parameter;
