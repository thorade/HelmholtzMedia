within HelmholtzMedia.Examples.ConvergenceTest;
model setSat
  package Medium = HelmholtzFluids.Carbondioxide_Short;
  Medium.SaturationProperties sat_T;
  Medium.SaturationProperties sat_p;
  Medium.SaturationProperties sat_dl;
  Medium.SaturationProperties sat_dv;
  Medium.DerPressureByTemperature dpT;
  Medium.DerTemperatureByPressure dTp;
  Medium.DerPressureByTemperature dpT2;
  Medium.DerTemperatureByPressure dTp2;

  Modelica.SIunits.TemperatureDifference T_err_abs;
  Modelica.SIunits.PressureDifference p_err_abs;
  Modelica.SIunits.Density dl_err_abs(min=-1e3) "Liquid density difference";
  Modelica.SIunits.Density dv_err_abs(min=-1e3) "Vapour density difference";

  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=8,
    startTime=1,
    height=Tcrit - Tmin +5,
    offset=Tmin)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  parameter Real eps = 1e-6;
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;

equation
  // forward
  sat_T = Medium.setSat_T(T=T_ramp.y);
  dpT = Medium.saturationPressure_derT(T=sat_T.Tsat);
  dTp = Medium.saturationTemperature_derp(p=sat_T.psat);
  dpT2 = 1/dTp;
  dTp2 = 1/dpT;

  // backward
  sat_p = Medium.setSat_p(p=sat_T.psat);
  sat_dl = Medium.setSat_d(d=sat_T.liq.d);
  sat_dv = Medium.setSat_d(d=sat_T.vap.d);

  T_err_abs = sat_T.Tsat - sat_p.Tsat;
  p_err_abs = sat_T.psat - sat_p.psat;
  dl_err_abs = sat_T.liq.d - sat_dl.liq.d;
  dv_err_abs = sat_T.vap.d - sat_dv.vap.d;

  annotation (experiment(
      StopTime=10,
      Tolerance=1e-005));
end setSat;
