within HelmholtzMedia.Examples.ConvergenceTest;
model TwoPhase_setState
  package Medium = HelmholtzFluids.Pentane;
  Medium.Density d;
  Medium.Temperature T;
  Medium.ThermodynamicState state;
  Medium.ThermodynamicState state_dT;
  Medium.ThermodynamicState state_pd;
  Medium.ThermodynamicState state_ph;
  Medium.ThermodynamicState state_ps;
  Medium.ThermodynamicState state_Ts;

Modelica.Blocks.Sources.Ramp T_sub(
    height=Tcrit - Tmin,
    offset=Tmin,
    duration=9,
    startTime=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  constant Medium.Density dmin=Medium.fluidLimits.DMIN;
  constant Medium.Density dcrit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Density dmax=Medium.fluidLimits.DMAX;
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;

equation
  d = dcrit;
  T = T_sub.y;

  // set state to valid two-phase values
  state=Medium.setState_Tx(T=T, x=0.5);

  // call setState functions
  state_dT=Medium.setState_dT(d=state.d, T=state.T, phase=0);
  state_pd=Medium.setState_pd(p=state.p, d=state.d, phase=0);
  state_ph=Medium.setState_ph(p=state.p, h=state.h, phase=0);
  state_ps=Medium.setState_ps(p=state.p, s=state.s, phase=0);
  state_Ts=Medium.setState_Ts(T=state.T, s=state.s, phase=0);

  annotation (experiment(StopTime=12, NumberOfIntervals=10000));
end TwoPhase_setState;
