within HelmholtzMedia.Examples.Parameter;
model State_dT_parameter "calculate state record from dT input"

  package Medium = HelmholtzFluids.R32;

  parameter Medium.Density d=1000;
  parameter Medium.Temperature T=304;

  Medium.ThermodynamicState state;
  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);

  // Medium.MassFraction x;
  // Medium.SurfaceTension sigma;
  // Medium.DynamicViscosity eta;
  // Medium.ThermalConductivity lambda;
  Medium.DerTemperatureByPressure dTp;
  Medium.DerPressureByTemperature dpT_sat;
  Medium.DerPressureByTemperature dpT;
  Medium.SpecificHeatCapacity cv;
  Medium.SpecificHeatCapacity cp;
  Medium.SpecificHeatCapacity cp0;
  Medium.IsentropicExponent gamma;
  Medium.VelocityOfSound a;
  Medium.DerTemperatureByPressure mu;
  Modelica.SIunits.IsothermalCompressibility kappa;
  Medium.IsobaricExpansionCoefficient beta;
  Medium.DerEnthalpyByPressure delta_T;

equation
  state=Medium.setState_dTX(d=d, T=T, phase=0);

  // x=Medium.vapourQuality(state);
  // sigma=Medium.surfaceTension(Medium.setSat_T(T=T));
  // eta=Medium.dynamicViscosity(state);
  // lambda=Medium.thermalConductivity(state);
  dTp=Medium.saturationTemperature_derp(p=state.p);
  dpT_sat=Medium.saturationPressure_derT(T=state.T);
  dpT=Medium.EoS.dpTd(f=f);
  cv=Medium.specificHeatCapacityCv(state=state);
  cp=Medium.specificHeatCapacityCp(state=state);
  cp0=f.R_s*(1-f.tau*f.tau*f.itt);
  gamma=Medium.isentropicExponent(state);
  a=Medium.velocityOfSound(state);
  mu=Medium.jouleThomsonCoefficient(state);
  kappa=Medium.isothermalCompressibility(state);
  beta=Medium.isobaricExpansionCoefficient(state);
  delta_T=Medium.isothermalThrottlingCoefficient(state);

  annotation (experiment(
      StopTime=2,
      Interval=1,
      Tolerance=1e-07));
end State_dT_parameter;
