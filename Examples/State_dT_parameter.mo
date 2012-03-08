within HelmholtzMedia.Examples;
model State_dT_parameter "calculate state record from dT input"

  package medium = HelmholtzFluids.R134a;

  parameter medium.Density d=700;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState state;

  // medium.MassFraction x;
  // medium.SurfaceTension sigma;
  // medium.DynamicViscosity eta;
  // medium.ThermalConductivity lambda;
  medium.SpecificHeatCapacity cv;
  Real dTp;
  Real dpT;

equation
  state=medium.setState_dTX(d=d, T=T, phase=0, X={1});

  // x=medium.vapourQuality(state);
  // sigma=medium.surfaceTension(medium.setSat_T(T=T));
  // eta=medium.dynamicViscosity(state);
  // lambda=medium.thermalConductivity(state);
  cv=medium.specificHeatCapacityCv(state=state);
  dpT=medium.saturationPressure_derT(T=state.T);
  dTp=medium.saturationTemperature_derp(p=state.p);

end State_dT_parameter;
