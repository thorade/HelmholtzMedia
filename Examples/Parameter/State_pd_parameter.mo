within HelmholtzMedia.Examples.Parameter;
model State_pd_parameter "calculate state record from pd input"

  package Medium = HelmholtzFluids.Propane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Density d=25;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_pd(p=p, d=d, phase=0);
end State_pd_parameter;
