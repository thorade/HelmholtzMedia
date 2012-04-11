within HelmholtzMedia.Examples;
model State_pT_parameter "calculate state record from pT input"

  package medium = HelmholtzMedia.HelmholtzFluids.Isopentane;

  parameter medium.AbsolutePressure p=101325;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState state;

  // pT always results in single phase states,
  // so this is a good place to test quantities that are defined for single phase only
  medium.Types.DerPressureByDensity dpdT;
  medium.Types.DerPressureByTemperature dpTd;
  medium.Types.DerEnthalpyByDensity dhdT;
  medium.Types.DerEnthalpyByTemperature dhTd;

  medium.DerDensityByTemperature ddTh;
  medium.DerDensityByPressure ddpT;
  medium.DerDensityByTemperature ddTp;
  medium.DerDensityByPressure ddph;
  medium.DerDensityByEnthalpy ddhp;

  medium.SpecificHeatCapacity cp;
  medium.SpecificHeatCapacity cv;
  medium.Types.DerTemperatureByPressure mu;
  medium.VelocityOfSound a;
  medium.ThermalConductivity lambda;
  medium.DynamicViscosity eta;
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
  mu=medium.jouleThomsonCoefficient(state);
  a=medium.velocityOfSound(state);

  lambda=medium.thermalConductivity(state);
  eta=medium.dynamicViscosity(state);
  Pr=medium.prandtlNumber(state);

end State_pT_parameter;
