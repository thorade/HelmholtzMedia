within HelmholtzMedia.Examples.ConvergenceTest;
model AncillaryFunctions
  package medium = HelmholtzFluids.Pentane;
  medium.Temperature Tsat;
  medium.AbsolutePressure psat;
  medium.Density dliq;
  medium.Density dvap;
  medium.SpecificEnthalpy hliq;
  medium.SpecificEntropy sliq;

  medium.Temperature T_p;
  medium.Temperature T_dl;
  medium.Temperature T_dv;
  medium.Temperature T_hl;
  medium.Temperature T_sl;

  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=5,
    height=T_crit - T_trip,
    offset=T_trip,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  constant medium.Temperature T_trip=medium.fluidConstants[1].triplePointTemperature;
  constant medium.Temperature T_crit=medium.fluidConstants[1].criticalTemperature;
  medium.EoS.HelmholtzDerivs fliq=medium.EoS.setHelmholtzDerivsFirst(d=dliq, T=Tsat);

equation
  // forward
  Tsat = T_ramp.y;
  psat = medium.Ancillary.saturationPressure_T(T=Tsat);
  dliq = medium.Ancillary.bubbleDensity_T(T=Tsat);
  dvap = medium.Ancillary.dewDensity_T(T=Tsat);
  hliq =  medium.EoS.h(f=fliq);
  sliq =  medium.EoS.s(f=fliq);

  // inverse
  T_p  = medium.Ancillary.saturationTemperature_p(p=psat);
  T_dl = medium.Ancillary.saturationTemperature_d(d=dliq);
  T_dv = medium.Ancillary.saturationTemperature_d(d=dvap);
  T_hl =  medium.Ancillary.saturationTemperature_h_liq(h=hliq);
  T_sl =  medium.Ancillary.saturationTemperature_s_liq(s=sliq);

annotation (experiment(
      StopTime=7,
      NumberOfIntervals=1000,
      Tolerance=1e-005));
end AncillaryFunctions;
