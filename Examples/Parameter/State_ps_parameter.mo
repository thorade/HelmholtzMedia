within HelmholtzMedia.Examples.Parameter;
model State_ps_parameter "calculate state record from ps input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=4.8813e+007;
  parameter Medium.SpecificEntropy s=593.986;

  Medium.ThermodynamicState inletState;

equation
  inletState=Interfaces.PartialHelmholtzMedium.setState_ps(p=p, s=s, phase=0);
end State_ps_parameter;
