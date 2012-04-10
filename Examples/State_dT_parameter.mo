within HelmholtzMedia.Examples;
model State_dT_parameter "calculate state record from dT input"

  package medium = HelmholtzFluids.Isopentane;

  parameter medium.Density d=33;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState state;

  // medium.MassFraction x;
  // medium.SurfaceTension sigma;
  // medium.DynamicViscosity eta;
  // medium.ThermalConductivity lambda;
  // medium.Types.DerTemperatureByPressure dTp;
  // medium.Types.DerPressureByTemperature dpT;
  medium.SpecificHeatCapacity cv;

equation
  state=medium.setState_dTX(d=d, T=T, phase=0, X={1});

  // x=medium.vapourQuality(state);
  // sigma=medium.surfaceTension(medium.setSat_T(T=T));
  // eta=medium.dynamicViscosity(state);
  // lambda=medium.thermalConductivity(state);
  // dpT=medium.saturationPressure_derT(T=state.T);
  // dTp=medium.saturationTemperature_derp(p=state.p);
  cv=medium.specificHeatCapacityCv(state=state);

end State_dT_parameter;
