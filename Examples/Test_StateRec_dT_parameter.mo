within HelmholtzFluids.Examples;
model Test_StateRec_dT_parameter
  package medium = HelmholtzFluids.Butane;
  parameter medium.Density d=228;
  parameter medium.Temperature T=298.15;
  medium.ThermodynamicState inlet;
  // medium.SpecificHeatCapacity cp;
  // medium.SpecificHeatCapacity cv;
  // medium.VelocityOfSound a;
  // medium.DynamicViscosity eta;
  // medium.ThermalConductivity lambda;
  // medium.PrandtlNumber Pr;
  medium.MassFraction x;
  medium.SurfaceTension sigma;

equation
  inlet=medium.setState_dTX(d=d, T=T, phase=0, X={1});

  // cp=medium.specificHeatCapacityCp(inlet);
  // cv=medium.specificHeatCapacityCv(inlet);
  // a=medium.velocityOfSound(inlet);

  // eta=medium.dynamicViscosity(inlet);
  // lambda=medium.thermalConductivity(inlet);

  // Pr=medium.prandtlNumber(inlet);

  x=medium.vapourQuality(inlet);
  sigma=medium.surfaceTension(medium.setSat_T(T=T));

end Test_StateRec_dT_parameter;
