within HelmholtzMedia.Examples;
model State_pT_parameter "calculate state record from pT input"

  package medium = HelmholtzFluids.Butane;

  parameter medium.AbsolutePressure p=101325;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState state;

  // pT always results in single phase states,
  // so this is a good place to test quantities that are defined for single phase only
  Interfaces.PartialHelmholtzMedium.Types.DerPressureByDensity pddT;
  Interfaces.PartialHelmholtzMedium.Types.DerPressureByTemperature pdTd;
  medium.SpecificHeatCapacity cp;
  medium.SpecificHeatCapacity cv;
  medium.VelocityOfSound a;
  medium.DynamicViscosity eta;
  medium.ThermalConductivity lambda;
  medium.PrandtlNumber Pr;

equation
  state=medium.setState_pTX(p=p, T=T, phase=0, X={1});

  pddT=medium.pressure_derd_T(state);
  pdTd=medium.pressure_derT_d(state);

  cp=medium.specificHeatCapacityCp(state);
  cv=medium.specificHeatCapacityCv(state);
  a=medium.velocityOfSound(state);

  eta=medium.dynamicViscosity(state);
  lambda=medium.thermalConductivity(state);

  Pr=medium.prandtlNumber(state);

end State_pT_parameter;
