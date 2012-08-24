within HelmholtzMedia.Examples.Sweep;
model setSat_sweep
  package Medium = HelmholtzFluids.Butane;
  Medium.SaturationProperties sat;

  Modelica.Blocks.Sources.Ramp Tramp(
    height=Tcrit - Tmin-(1e-12),
    offset=Tmin,
    duration=8,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Ramp pramp(
    duration=8,
    startTime=1,
    height=pcrit - pmin - 1,
    offset=pmin + 1)
                 annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

equation
// sat = Medium.setSat_T(T=Tramp.y);
sat = Medium.setSat_p(p=pramp.y);

end setSat_sweep;
