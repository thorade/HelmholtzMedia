within HelmholtzMedia.Examples.Parameter;
model State_Ts_parameter "calculate state record from Ts input"

  package Medium = HelmholtzFluids.Propane;

  parameter Medium.Temperature T=300.05;
  parameter Medium.SpecificEntropy s=2353.05;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_Ts(T=T, s=s, phase=0);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end State_Ts_parameter;
