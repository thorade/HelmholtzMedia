within HelmholtzMedia.Examples.Parameter;
model State_dT_parameter "calculate state record from dT input"

  package Medium = HelmholtzFluids.HMDS;

  parameter Medium.Density d=1e-3;
  parameter Medium.Temperature T=298.15;

  Medium.ThermodynamicState state;
  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);

  // Medium.MassFraction x;
  // Medium.SurfaceTension sigma;
  // Medium.DynamicViscosity eta;
  // Medium.ThermalConductivity lambda;
  // Medium.Types.DerTemperatureByPressure dTp;
  // Medium.Types.DerPressureByTemperature dpT;
  Medium.SpecificHeatCapacity cv;
  Medium.SpecificHeatCapacity cp;
  Medium.SpecificHeatCapacity cp0;
  Medium.IsentropicExponent gamma;
  Medium.VelocityOfSound a;
  Medium.Types.DerTemperatureByPressure mu;
  Modelica.SIunits.IsothermalCompressibility kappa;
  Medium.IsobaricExpansionCoefficient beta;
  Medium.DerEnthalpyByPressure delta_T;

equation
  state=Interfaces.PartialHelmholtzMedium.setState_dT(d=d, T=T, phase=0);

  // x=Medium.vapourQuality(state);
  // sigma=Medium.surfaceTension(Medium.setSat_T(T=T));
  // eta=Medium.dynamicViscosity(state);
  // lambda=Medium.thermalConductivity(state);
  // dpT=Medium.saturationPressure_derT(T=state.T);
  // dTp=Medium.saturationTemperature_derp(p=state.p);
  cv=Medium.specificHeatCapacityCv(state=state);
  cp=Medium.specificHeatCapacityCp(state=state);
  cp0=f.R*(1-f.tau*f.tau*f.itt);
  gamma=Medium.isentropicExponent(state);
  a=Medium.velocityOfSound(state);
  mu=Medium.jouleThomsonCoefficient(state);
  kappa=Medium.isothermalCompressibility(state);
  beta=Medium.isobaricExpansionCoefficient(state);
  delta_T=Medium.isothermalThrottlingCoefficient(state);

end State_dT_parameter;
