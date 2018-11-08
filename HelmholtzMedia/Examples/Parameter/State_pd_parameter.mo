within HelmholtzMedia.Examples.Parameter;
model State_pd_parameter "calculate state record from pd input"

  package Medium = HelmholtzFluids.Propane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Density d=1e-6;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_pd(p=p, d=d, phase=0);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end State_pd_parameter;
