within HelmholtzMedia.Examples.Sweep;
model setSat_sweep
  package medium = HelmholtzFluids.Butane;
  medium.SaturationProperties sat;
  medium.Temperature T_test;

  Modelica.Blocks.Sources.Ramp Tramp(
    duration=8,
    startTime=1,
    height=Tcrit - Tmin - 1,
    offset=Tmin + 1)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Ramp pramp(
    duration=8,
    startTime=1,
    height=pcrit - pmin - 1,
    offset=pmin + 1)
                 annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

protected
  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature Tmax=medium.fluidLimits.TMAX;
  constant medium.AbsolutePressure pmin=medium.fluidLimits.PMIN;
  constant medium.AbsolutePressure pcrit=medium.fluidConstants[1].criticalPressure;
  constant medium.AbsolutePressure pmax=medium.fluidLimits.PMAX;

equation
sat = medium.setSat_T(T=Tramp.y);
// sat = medium.setSat_p(p=pramp.y);

T_test = Tramp.y - medium.Ancillary.saturationTemperature_d(sat.vap.d);
  annotation (experiment(StopTime=11), __Dymola_experimentSetupOutput);
end setSat_sweep;
