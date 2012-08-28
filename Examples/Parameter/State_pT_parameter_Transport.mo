within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter_Transport "calculate state record from pT input"

  package Medium = HelmholtzFluids.Propane;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;

  // pT always results in single phase states,
  // so this is a good place to test quantities that are defined for single phase only
  Medium.Types.DerPressureByDensity dpdT;
  Medium.Types.DerPressureByTemperature dpTd;
  Medium.Types.DerEnthalpyByDensity dhdT;
  Medium.Types.DerEnthalpyByTemperature dhTd;

  Medium.DerDensityByTemperature ddTh;
  Medium.DerDensityByPressure ddpT;
  Medium.DerDensityByTemperature ddTp;
  Medium.DerDensityByPressure ddph;
  Medium.DerDensityByEnthalpy ddhp;

  Medium.SpecificHeatCapacity cp;
  Medium.SpecificHeatCapacity cv;
  Medium.Types.DerTemperatureByPressure mu;
  Medium.VelocityOfSound a;
  Medium.ThermalConductivity lambda;
  Medium.DynamicViscosity eta;
  Medium.PrandtlNumber Pr;

equation
  state=Medium.setState_pTX(p=p, T=T, phase=0, X={1});

  dpdT=Medium.pressure_derd_T(state);
  dpTd=Medium.pressure_derT_d(state);
  dhdT=Medium.specificEnthalpy_derd_T(state);
  dhTd=Medium.specificEnthalpy_derT_d(state);

  ddTh=Medium.density_derT_h(state);
  ddpT=Medium.density_derp_T(state);
  ddTp=Medium.density_derT_p(state);
  ddph=Medium.density_derp_h(state);
  ddhp=Medium.density_derh_p(state);

  cp=Medium.specificHeatCapacityCp(state);
  cv=Medium.specificHeatCapacityCv(state);
  mu=Medium.jouleThomsonCoefficient(state);
  a=Medium.velocityOfSound(state);

  lambda=Medium.thermalConductivity(state);
  eta=Medium.dynamicViscosity(state);
  Pr=Medium.prandtlNumber(state);

end State_pT_parameter_Transport;
