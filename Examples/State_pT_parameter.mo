within HelmholtzFluids.Examples;
model State_pT_parameter "calculate state record from pT input"

  package medium = HelmholtzFluids.Butane;

  parameter medium.AbsolutePressure p=101325;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState inletState;

  // pT always results in single phase states,
  // so this is a good place to test quantities that are defined for single phase only
  medium.SpecificHeatCapacity cp;
  medium.SpecificHeatCapacity cv;
  medium.VelocityOfSound a;
  medium.DynamicViscosity eta;
  medium.ThermalConductivity lambda;
  medium.PrandtlNumber Pr;

equation
  inletState=medium.setState_pTX(p=p, T=T, phase=0, X={1});

  cp=medium.specificHeatCapacityCp(inletState);
  cv=medium.specificHeatCapacityCv(inletState);
  a=medium.velocityOfSound(inletState);

  eta=medium.dynamicViscosity(inletState);
  lambda=medium.thermalConductivity(inletState);

  Pr=medium.prandtlNumber(inletState);

end State_pT_parameter;
