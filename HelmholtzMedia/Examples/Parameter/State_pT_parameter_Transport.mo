within HelmholtzMedia.Examples.Parameter;
model State_pT_parameter_Transport "calculate state record from pT input"

  package Medium = HelmholtzFluids.Carbondioxide;

  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;
  Medium.ThermodynamicState state;
  // pT always results in single phase states

  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);

  // order of properties as in RefProp
  // derived properties
  Medium.SpecificHeatCapacity cv;
  Medium.SpecificHeatCapacity cp;
  Medium.SpecificHeatCapacity cp0;
  Medium.IsentropicExponent gamma;
  Medium.VelocityOfSound a;
  Medium.DerTemperatureByPressure mu;
  // transport proerties
  Medium.ThermalConductivity lambda;
  Medium.DynamicViscosity eta;
  //Medium.PrandtlNumber Pr;
  // more derived properties
  Modelica.SIunits.IsothermalCompressibility kappa;
  Medium.IsobaricExpansionCoefficient beta;
  Medium.DerEnthalpyByPressure delta_T;

equation
  state=Medium.setState_pTX(p=p, T=T, phase=0, X={1});

  // derived properties
  cv=Medium.specificHeatCapacityCv(state);
  cp=Medium.specificHeatCapacityCp(state);
  cp0=f.R*(1-f.tau*f.tau*f.itt);
  gamma=Medium.isentropicExponent(state);
  a=Medium.velocityOfSound(state);
  mu=Medium.jouleThomsonCoefficient(state);
  // transport properties
  lambda=Medium.thermalConductivity(state);
  eta=Medium.dynamicViscosity(state);
  //Pr=Medium.prandtlNumber(state);
  // more derived properties
  kappa=Medium.isothermalCompressibility(state);
  beta=Medium.isobaricExpansionCoefficient(state);
  delta_T=Medium.isothermalThrottlingCoefficient(state);

  annotation (experiment(Tolerance=1e-06));
end State_pT_parameter_Transport;
