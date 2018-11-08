within HelmholtzMedia.Examples.Parameter;
model State_ph_parameter "calculate state record from ph input"

  package Medium = HelmholtzFluids.Helium;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.SpecificEnthalpy h=2e3;

  Medium.ThermodynamicState inletState;

equation
  inletState=Medium.setState_phX(p=p, h=h, phase=0);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end State_ph_parameter;
