within HelmholtzMedia.Examples;
model State_pT_parameter2 "calculate state record from pT input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=3796000;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;
  // Medium.SpecificEnthalpy h;

equation
  state=Medium.setState_pTX(p=p, T=T, phase=0);
  // h = Medium.specificEnthalpy_pT(p=p, T=T);
end State_pT_parameter2;
