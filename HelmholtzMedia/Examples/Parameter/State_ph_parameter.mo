within HelmholtzMedia.Examples.Parameter;
model State_ph_parameter "calculate state record from ph input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.SpecificEnthalpy h=2e3;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_phX(p=p, h=h, phase=0);
end State_ph_parameter;
