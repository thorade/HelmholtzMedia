within HelmholtzFluids.Examples;
model State_dT_parameter "calculate state record from dT input"

  package medium = HelmholtzFluids.Butane;

  parameter medium.Density d=228;
  parameter medium.Temperature T=300;

  medium.ThermodynamicState inletState;

  // medium.MassFraction x;
  // medium.SurfaceTension sigma;
  // medium.DynamicViscosity eta;
  // medium.ThermalConductivity lambda;
  medium.SpecificHeatCapacity cv;

equation
  inletState=medium.setState_dTX(d=d, T=T, phase=0, X={1});

  // x=medium.vapourQuality(inletState);
  // sigma=medium.surfaceTension(medium.setSat_T(T=T));
  // eta=medium.dynamicViscosity(inletState);
  // lambda=medium.thermalConductivity(inletState);
  cv=medium.specificHeatCapacityCv(inletState);

end State_dT_parameter;
