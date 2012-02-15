within HelmholtzFluids.Examples;
model State_ph_parameter "calculate state record from ph input"
  package medium = HelmholtzFluids.Butane;

  parameter medium.AbsolutePressure p=101325;
  parameter medium.SpecificEnthalpy h=10e3;

  medium.ThermodynamicState inletState;

equation
  inletState=medium.setState_phX(p=p, h=h, phase=0, X={1});
end State_ph_parameter;
