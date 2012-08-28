within HelmholtzMedia.Examples.Parameter;
model setSat_parameter
  package medium = HelmholtzFluids.Propane;

  medium.SaturationProperties sat;
  medium.SaturationProperties sat_p;
  medium.SaturationProperties sat_dl;
  medium.SaturationProperties sat_dv;

  parameter medium.Temperature T_input= 298.15;

protected
  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;

equation
  // forward
  sat = medium.setSat_T(T=T_input);

  // backward
  sat_p  = medium.setSat_p(p=sat.psat);
  sat_dl = medium.setSat_d(d=sat.liq.d);
  sat_dv = medium.setSat_d(d=sat.vap.d);

  annotation (experiment(NumberOfIntervals=1));
end setSat_parameter;
