within HelmholtzFluids.Examples;
model Test_StateRec_ph_parameter
  package medium = HelmholtzFluids.Butane;
  parameter medium.AbsolutePressure p=2e6;
  parameter medium.SpecificEnthalpy h=10e3;
  medium.ThermodynamicState inletState;
  medium.Density d_in;
  medium.Temperature T_in;
  medium.SaturationProperties test;

equation
  inletState=medium.setState_phX(p=p, h=h, phase=0, X={1});
  test= medium.setSat_p(p=p);
  d_in=inletState.d;
  T_in=inletState.T;
end Test_StateRec_ph_parameter;
