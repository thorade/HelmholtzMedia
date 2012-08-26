within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter2 "calculate state record from pT input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=2.6e+006;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;
  // Medium.SpecificEnthalpy h;

equation
  state=Medium.setState_pTX(p=p, T=T, phase=0);
  // h = Medium.specificEnthalpy_pT(p=p, T=T);
end State_pT_parameter2;
