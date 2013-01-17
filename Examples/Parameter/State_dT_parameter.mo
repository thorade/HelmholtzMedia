within HelmholtzMedia.Examples.Parameter;
model State_dT_parameter "calculate state record from dT input"

  package Medium = HelmholtzFluids.R134a;

  parameter Medium.Density d=1e-3;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;

  // Medium.MassFraction x;
  // Medium.SurfaceTension sigma;
  // Medium.DynamicViscosity eta;
  Medium.ThermalConductivity lambda;
  // Medium.Types.DerTemperatureByPressure dTp;
  // Medium.Types.DerPressureByTemperature dpT;
  Medium.SpecificHeatCapacity cv;
  Medium.SpecificHeatCapacity cp;
  // Medium.VelocityOfSound a;

equation
  state=Medium.setState_dTX(d=d, T=T, phase=0);

  // x=Medium.vapourQuality(state);
  // sigma=Medium.surfaceTension(Medium.setSat_T(T=T));
  // eta=Medium.dynamicViscosity(state);
  lambda=Medium.thermalConductivity(state);
  // dpT=Medium.saturationPressure_derT(T=state.T);
  // dTp=Medium.saturationTemperature_derp(p=state.p);
  cv=Medium.specificHeatCapacityCv(state=state);
  cp=Medium.specificHeatCapacityCp(state=state);
  // a=Medium.velocityOfSound(state=state);

end State_dT_parameter;
