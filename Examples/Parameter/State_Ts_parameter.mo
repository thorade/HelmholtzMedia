within HelmholtzMedia.Examples.Parameter;
model State_Ts_parameter "calculate state record from Ts input"

  package Medium = HelmholtzFluids.R134a;

  parameter Medium.Temperature T=379.747;
  parameter Medium.SpecificEntropy s=2e3;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_Ts(T=T, s=s, phase=0);
end State_Ts_parameter;
