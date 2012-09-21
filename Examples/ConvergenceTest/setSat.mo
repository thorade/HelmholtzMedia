within HelmholtzMedia.Examples.ConvergenceTest;
model setSat
  package medium = HelmholtzFluids.Pentane;
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

equation
  // forward
  sat = medium.setSat_T(T=T_ramp.y);

  // backward
  sat_p = medium.setSat_p(p=sat.psat);
  sat_dl = medium.setSat_d(d=sat.liq.d);
  sat_dv = medium.setSat_d(d=sat.vap.d);

  annotation (experiment(
      StopTime=10,
      NumberOfIntervals=10000,
      Tolerance=1e-005));
end setSat;
