within HelmholtzMedia.Examples.ConvergenceTest;
model Ancillary_StartValues
  package Medium = HelmholtzFluids.Carbondioxide_Short;
  Medium.AbsolutePressure p;
  Medium.Temperature T;
  Medium.Density d;
  // ========
  Medium.AbsolutePressure psat;
  // ========
  Medium.AbsolutePressure p_dT_Waals;
  Medium.Density d_pT_Soave;
  Medium.Temperature T_pd_Waals;

Modelica.Blocks.Sources.Ramp p_sub(
    duration=4,
    startTime=0.1,
    height=pcrit - pmin,
    offset=pmin)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

Modelica.Blocks.Sources.Ramp p_super(
    duration=5,
    startTime=6,
    offset=0,
    height=pmax - pcrit)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

  Modelica.Blocks.Sources.Sine T_sine(
    startTime=0,
    amplitude=(Tmax - Tmin)/2,
    offset=(Tmax - Tmin)/2 + Tmin,
    f=1) annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure pmin=1e-6;//Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

equation
  p = p_sub.y + p_super.y;
  T = T_sine.y;

  // get d and psat from EoS
  d=Medium.density_pT(p=p, T=T);
  psat=Medium.saturationPressure(T);

  // call Ancillary start value functions
  d_pT_Soave = Medium.Ancillary.density_pT_Soave(T=T, p=p, psat=psat);
  p_dT_Waals = Medium.Ancillary.pressure_dT_Waals(d=d, T=T);
  T_pd_Waals = Medium.Ancillary.temperature_pd_Waals(p=p, d=d);

  annotation (experiment(StopTime=12));
end Ancillary_StartValues;
