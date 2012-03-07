within HelmholtzMedia.Examples;
model State_ps_parameter "calculate state record from ps input"

  package medium = HelmholtzFluids.Butane;

  parameter medium.AbsolutePressure p=101325;
  parameter medium.SpecificEntropy s=1e3;

  medium.ThermodynamicState inletState;

equation
  inletState=medium.setState_psX(p=p, s=s, phase=0, X={1});
end State_ps_parameter;
