within HelmholtzMedia.Validation.Convergence;
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
  Medium.Internal.EoS.HelmholtzDerivs fliq=Medium.Internal.EoS.setHelmholtzDerivs(d=dliq, T=Tsat, order=1);

equation
  // forward
  Tsat = T_ramp.y;
  psat = Medium.Internal.Ancillary.saturationPressure_T(T=Tsat);
  pmelt = Medium.Internal.Ancillary.meltingPressure_T(T=Tsat);
  dliq = Medium.Internal.Ancillary.bubbleDensity_T(T=Tsat);
  dvap = Medium.Internal.Ancillary.dewDensity_T(T=Tsat);
  hliq =  Medium.Internal.EoS.h(f=fliq);
  sliq =  Medium.Internal.EoS.s(f=fliq);

  // inverse
  T_ps = Medium.Internal.Ancillary.saturationTemperature_p(p=psat);
  T_pm = Medium.Internal.Ancillary.meltingTemperature_p(p=pmelt);
  T_dl = Medium.Internal.Ancillary.saturationTemperature_d(d=dliq);
  T_dv = Medium.Internal.Ancillary.saturationTemperature_d(d=dvap);
  T_hl =  Medium.Internal.Ancillary.saturationTemperature_h_liq(h=hliq);
  T_sl =  Medium.Internal.Ancillary.saturationTemperature_s_liq(s=sliq);

annotation (experiment(
      StopTime=7,
      Tolerance=1e-005));
end Ancillary_Saturation;
