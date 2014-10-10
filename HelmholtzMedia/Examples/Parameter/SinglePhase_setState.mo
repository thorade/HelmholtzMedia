within HelmholtzMedia.Examples.Parameter;
model SinglePhase_setState
  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;
  Medium.ThermodynamicState state_dT;
  Medium.ThermodynamicState state_pd;
  Medium.ThermodynamicState state_ph;
  Medium.ThermodynamicState state_ps;
  Medium.ThermodynamicState state_Ts;

equation
  // set state to valid single-phase values
  state=Medium.setState_pT(p=p, T=T, phase=0);

  // call other setState functions
  state_dT=Medium.setState_dT(d=state.d, T=state.T, phase=0);
  state_pd=Medium.setState_pd(p=state.p, d=state.d, phase=0);
  state_ph=Medium.setState_ph(p=state.p, h=state.h, phase=0);
  state_ps=Medium.setState_ps(p=state.p, s=state.s, phase=0);
  state_Ts=Medium.setState_Ts(T=state.T, s=state.s, phase=0);

  annotation (experiment(StopTime=12));
end SinglePhase_setState;
