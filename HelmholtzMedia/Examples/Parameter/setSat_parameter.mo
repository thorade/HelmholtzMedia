within HelmholtzMedia.Examples.Parameter;
model setSat_parameter
  package Medium = HelmholtzFluids.Helium;

  Medium.SaturationProperties sat;
  Medium.SaturationProperties sat_p;
  Medium.SaturationProperties sat_dl;
  Medium.SaturationProperties sat_dv;

  parameter Medium.Temperature T_input= Tcrit - 0.1;

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;

equation
  // forward
  sat = Medium.setSat_T(T=T_input);

  // backward
  sat_p  = Medium.setSat_p(p=sat.psat);
  sat_dl = Medium.setSat_d(d=sat.liq.d);
  sat_dv = Medium.setSat_d(d=sat.vap.d);

end setSat_parameter;
