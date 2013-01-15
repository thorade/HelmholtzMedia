within HelmholtzMedia.Examples.Validation;
model MaxwellLoop "show Maxwell Loops"
  package Medium = HelmholtzFluids.Butane;
  parameter Medium.Temperature T = 0.7*Tcrit;

  Medium.Density d;
  Medium.ThermodynamicState state = Medium.setState_dTX(d=d, T=T);
  Medium.EoS.HelmholtzDerivs f = Medium.EoS.setHelmholtzDerivsThird(d=d, T=T, phase=1);
  Medium.AbsolutePressure p = Medium.EoS.p(f);

protected
  constant Medium.MolarMass MM = Medium.fluidConstants[1].molarMass;
  constant Medium.SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";

  constant Medium.Density dcrit=MM/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Density dmin=Medium.fluidLimits.DMIN;
  constant Medium.Density dmax=Medium.fluidLimits.DMAX;

  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;

public
Modelica.Blocks.Sources.Ramp Ramp_dvap(
    duration=4,
    startTime=0.1,
    height=dcrit - dmin,
    offset=dmin)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
Modelica.Blocks.Sources.Ramp Ramp_dliq(
    duration=5,
    startTime=6,
    height=dmax - dcrit,
    offset=0)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

equation
  d = max({Ramp_dvap.y, 0.9*Medium.Ancillary.dewDensity_T(T=T)}) + min({Ramp_dliq.y, 1.05*Medium.Ancillary.bubbleDensity_T(T=T)});

  annotation (experiment(StopTime=12, __Dymola_NumberOfIntervals=1000),
                                      __Dymola_experimentSetupOutput);
end MaxwellLoop;
