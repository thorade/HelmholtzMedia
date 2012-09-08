within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter_Transport "calculate state record from pT input"

  package Medium = HelmholtzFluids.Propane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;
  Medium.ThermodynamicState state; // pT always results in single phase states

  // derived properties
  Medium.SpecificHeatCapacity cp;
  Medium.SpecificHeatCapacity cv;
  Medium.Types.DerTemperatureByPressure mu;
  Medium.VelocityOfSound a;

  // transport proerties and derived properties
  Medium.ThermalConductivity lambda;
  Medium.DynamicViscosity eta;
  Medium.PrandtlNumber Pr;

equation
  state=Medium.setState_pTX(p=p, T=T, phase=0, X={1});

  cp=Medium.specificHeatCapacityCp(state);
  cv=Medium.specificHeatCapacityCv(state);
  mu=Medium.jouleThomsonCoefficient(state);
  a=Medium.velocityOfSound(state);

  lambda=Medium.thermalConductivity(state);
  eta=Medium.dynamicViscosity(state);
  Pr=Medium.prandtlNumber(state);

end State_pT_parameter_Transport;
