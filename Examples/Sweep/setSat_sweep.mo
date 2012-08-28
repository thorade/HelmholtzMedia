within HelmholtzMedia.Examples.Sweep;
model setSat_sweep
  package medium = HelmholtzFluids.Butane;
  medium.SaturationProperties sat;
  medium.SaturationProperties sat_p;
  medium.SaturationProperties sat_dl;
  medium.SaturationProperties sat_dv;

  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=8,
    startTime=1,
    height=Tcrit - Tmin,
    offset=Tmin)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature Tmax=medium.fluidLimits.TMAX;

equation
  sat = medium.setSat_T(T=T_ramp.y);
  sat_p = medium.setSat_p(p=sat.psat);
  sat_dl = medium.setSat_d(d=sat.liq.d);
  sat_dv = medium.setSat_d(d=sat.vap.d);

  annotation (experiment(
      StopTime=11,
      NumberOfIntervals=1000,
      Tolerance=1e-005),               __Dymola_experimentSetupOutput);
end setSat_sweep;
