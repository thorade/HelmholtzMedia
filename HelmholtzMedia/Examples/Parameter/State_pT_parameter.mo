within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter "calculate state record from pT input"

  package Medium = HelmholtzFluids.Helium;

  parameter Medium.AbsolutePressure p=986636; //101325;
  parameter Medium.Temperature T=3.72952; //2.1768;

  Medium.ThermodynamicState state;
  //Medium.ThermodynamicState state_ph;
  Medium.ThermodynamicState state_ps;

equation
  state=Medium.setState_pT(p=p, T=T, phase=0);
  //state_ph=Medium.setState_ph(p=Medium.pressure(state), h=Medium.specificEnthalpy(state), phase=0);
  state_ps=Medium.setState_ps(p=Medium.pressure(state), s=Medium.specificEntropy(state), phase=0);
end State_pT_parameter;
