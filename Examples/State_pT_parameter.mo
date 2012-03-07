within HelmholtzMedia.Examples;
model State_pT_parameter "calculate state record from pT input"

  package medium = HelmholtzFluids.Butane;

  parameter medium.AbsolutePressure p=101325;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState state;

  // pT always results in single phase states,
  // so this is a good place to test quantities that are defined for single phase only
  Interfaces.PartialHelmholtzMedium.Types.DerPressureByDensity dpdT;
  Interfaces.PartialHelmholtzMedium.Types.DerPressureByTemperature dpTd;
  Interfaces.PartialHelmholtzMedium.Types.DerEnthalpyByDensity dhdT;
  Interfaces.PartialHelmholtzMedium.Types.DerEnthalpyByTemperature dhTd;

  Interfaces.PartialHelmholtzMedium.DerDensityByTemperature ddTh;
  Interfaces.PartialHelmholtzMedium.DerDensityByPressure ddpT;
  Interfaces.PartialHelmholtzMedium.DerDensityByTemperature ddTp;
  Interfaces.PartialHelmholtzMedium.DerDensityByPressure ddph;
  Interfaces.PartialHelmholtzMedium.DerDensityByEnthalpy ddhp;

  medium.SpecificHeatCapacity cp;
  medium.SpecificHeatCapacity cv;
  medium.VelocityOfSound a;
  medium.DynamicViscosity eta;
  medium.ThermalConductivity lambda;
  medium.PrandtlNumber Pr;

equation
  state=medium.setState_pTX(p=p, T=T, phase=0, X={1});

  dpdT=medium.pressure_derd_T(state);
  dpTd=medium.pressure_derT_d(state);
  dhdT=medium.specificEnthalpy_derd_T(state);
  dhTd=medium.specificEnthalpy_derT_d(state);

  ddTh=medium.density_derT_h(state);
  ddpT=medium.density_derp_T(state);
  ddTp=medium.density_derT_p(state);
  ddph=medium.density_derp_h(state);
  ddhp=medium.density_derh_p(state);

  cp=medium.specificHeatCapacityCp(state);
  cv=medium.specificHeatCapacityCv(state);
  a=medium.velocityOfSound(state);

  eta=medium.dynamicViscosity(state);
  lambda=medium.thermalConductivity(state);

  Pr=medium.prandtlNumber(state);

end State_pT_parameter;
