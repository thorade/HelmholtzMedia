within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter "calculate state record from pT input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=3.79723e+006;
  parameter Medium.Temperature T=425.125;

  Medium.ThermodynamicState state;

equation
  state=Medium.setState_pT(p=p, T=T, phase=0);
end State_pT_parameter;
