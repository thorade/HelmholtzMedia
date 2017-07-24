within HelmholtzMedia.Examples.ConvergenceTest;
model Ancillary_Saturation
  package Medium = HelmholtzFluids.Carbondioxide_Short;
  Medium.Temperature Tsat;
  Medium.AbsolutePressure psat;
  Medium.AbsolutePressure pmelt;
  Medium.Density dliq;
  Medium.Density dvap;
  Medium.SpecificEnthalpy hliq;
  Medium.SpecificEntropy sliq;

  Medium.Temperature T_ps;
  Medium.Temperature T_pm;
  Medium.Temperature T_dl;
  Medium.Temperature T_dv;
  Medium.Temperature T_hl;
  Medium.Temperature T_sl;

  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=5,
    height=T_crit - T_trip,
    offset=T_trip,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  constant Medium.Temperature T_trip=Medium.fluidConstants[1].triplePointTemperature;
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  Medium.EoS.HelmholtzDerivs fliq=Medium.EoS.setHelmholtzDerivsFirst(d=dliq, T=Tsat);

equation
  // forward
  Tsat = T_ramp.y;
  psat = Medium.Ancillary.saturationPressure_T(T=Tsat);
  pmelt = Medium.Ancillary.meltingPressure_T(T=Tsat);
  dliq = Medium.Ancillary.bubbleDensity_T(T=Tsat);
  dvap = Medium.Ancillary.dewDensity_T(T=Tsat);
  hliq =  Medium.EoS.h(f=fliq);
  sliq =  Medium.EoS.s(f=fliq);

  // inverse
  T_ps = Medium.Ancillary.saturationTemperature_p(p=psat);
  T_pm = Medium.Ancillary.meltingTemperature_p(p=pmelt);
  T_dl = Medium.Ancillary.saturationTemperature_d(d=dliq);
  T_dv = Medium.Ancillary.saturationTemperature_d(d=dvap);
  T_hl =  Medium.Ancillary.saturationTemperature_h_liq(h=hliq);
  T_sl =  Medium.Ancillary.saturationTemperature_s_liq(s=sliq);

annotation (experiment(
      StopTime=7,
      Tolerance=1e-005));
end Ancillary_Saturation;
