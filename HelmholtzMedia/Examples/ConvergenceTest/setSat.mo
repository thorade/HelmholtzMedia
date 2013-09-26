within HelmholtzMedia.Examples.ConvergenceTest;
model setSat
  package Medium = HelmholtzFluids.Helium;
  Medium.SaturationProperties sat;
  Medium.SaturationProperties sat_p;
  Medium.SaturationProperties sat_dl;
  Medium.SaturationProperties sat_dv;
  Medium.Types.DerPressureByTemperature dpT;
  Medium.Types.DerTemperatureByPressure dTp;

  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=8,
    startTime=1,
    height=Tcrit - Tmin +5,
    offset=Tmin)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;

equation
  // forward
  sat = Medium.setSat_T(T=T_ramp.y);
  dpT = Medium.saturationPressure_derT(T=sat.Tsat);
  dTp = Medium.saturationTemperature_derp(p=sat.psat);

  // backward
  sat_p = Medium.setSat_p(p=sat.psat);
  sat_dl = Medium.setSat_d(d=sat.liq.d);
  sat_dv = Medium.setSat_d(d=sat.vap.d);

  annotation (experiment(
      StopTime=10,
      NumberOfIntervals=10000,
      Tolerance=1e-005));
end setSat;
