within HelmholtzMedia.Examples.Sweep;
model setSat_sweep
  package medium = HelmholtzFluids.Butane;
  medium.SaturationProperties sat;

  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=8,
    startTime=1,
    height=Tcrit - Tmin - 1,
    offset=Tmin + 1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Ramp p_ramp(
    duration=8,
    startTime=1,
    height=pcrit - pmin - 1,
    offset=pmin + 1)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Ramp d_vap_ramp(
    duration=8,
    startTime=1,
    height=dcrit - dmin - 1e-3,
    offset=dmin + 1e-3)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Modelica.Blocks.Sources.Ramp d_liq_ramp(
    duration=8,
    startTime=1,
    offset=dcrit,
    height=dmax - dcrit - 100)
    annotation (Placement(transformation(extent={{-80,-80},{-60,-60}})));

protected
  constant medium.MolarMass MM = medium.fluidConstants[1].molarMass;
  constant medium.SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant medium.Density dcrit=MM/medium.fluidConstants[1].criticalMolarVolume;
  constant medium.Density dmin=medium.fluidLimits.DMIN;
  constant medium.Density dmax=medium.fluidLimits.DMAX;

  constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
  constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature Tmax=medium.fluidLimits.TMAX;

  constant medium.AbsolutePressure pmin=medium.fluidLimits.PMIN;
  constant medium.AbsolutePressure pcrit=medium.fluidConstants[1].criticalPressure;
  constant medium.AbsolutePressure pmax=medium.fluidLimits.PMAX;

equation
// sat = medium.setSat_T(T=T_ramp.y);
// sat = medium.setSat_p(p=p_ramp.y);
sat = medium.setSat_d(d=d_vap_ramp.y);

  annotation (experiment(StopTime=11), __Dymola_experimentSetupOutput);
end setSat_sweep;
