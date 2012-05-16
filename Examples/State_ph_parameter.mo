within HelmholtzMedia.Examples;
model State_ph_parameter "calculate state record from ph input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=16e5;
  parameter Medium.SpecificEnthalpy h=500e3;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_phX(p=p, h=h, phase=0, X=fill(0,0));
end State_ph_parameter;
