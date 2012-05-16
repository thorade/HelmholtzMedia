within HelmholtzMedia.Examples;
model setSat_sweep
  package Medium = HelmholtzFluids.R134a;
  Medium.SaturationProperties sat;
//Medium.ThermodynamicState liq;
//Medium.ThermodynamicState vap;

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

  Modelica.Blocks.Sources.Ramp Tramp(
    height=Tcrit - Tmin-(1e-12),
    offset=Tmin,
    duration=8,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.Ramp pramp(
    height=pcrit - pmin,
    offset=pmin,
    duration=8,
    startTime=1) annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

equation
  sat = Medium.setSat_T(T=Tramp.y);
//sat = Medium.setSat_p(p=pramp.y);

//liq = sat.liq;
//vap = sat.vap;

  annotation (experiment(
      StopTime=10,
      NumberOfIntervals=10000,
      Tolerance=1e-008),               __Dymola_experimentSetupOutput);
end setSat_sweep;
