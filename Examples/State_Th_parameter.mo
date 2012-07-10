within HelmholtzMedia.Examples;
model State_Th_parameter "calculate state record from Th input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.Temperature T=298.15;
  parameter Medium.SpecificEnthalpy h=500e3;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_ThX(T=T, h=h, phase=0);
end State_Th_parameter;
