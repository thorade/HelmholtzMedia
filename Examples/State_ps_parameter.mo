within HelmholtzMedia.Examples;
model State_ps_parameter "calculate state record from ps input"

  package Medium = HelmholtzFluids.Butane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.SpecificEntropy s=1e3;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_psX(p=p, s=s, phase=0);
end State_ps_parameter;
