within HelmholtzFluids.Examples;
model Test_StateRec_pT_parameter
  package medium = HelmholtzFluids.R134a;
  parameter medium.AbsolutePressure p=101325;
  parameter medium.Temperature T=298.15;
  medium.ThermodynamicState inletState;

equation
  inletState=medium.setState_pTX(p=p, T=T, phase=0, X={1});

end Test_StateRec_pT_parameter;
