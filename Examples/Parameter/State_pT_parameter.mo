within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter "calculate state record from pT input"

  package Medium = HelmholtzFluids.Propane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;

equation
  state=Medium.setState_pT(p=p, T=T, phase=0);
end State_pT_parameter;
