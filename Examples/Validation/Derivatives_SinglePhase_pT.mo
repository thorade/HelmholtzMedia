within HelmholtzMedia.Examples.Validation;
model Derivatives_SinglePhase_pT
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Butane;
  // p and T always result in single-phase
  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;

// Density derivatives
  Medium.DerDensityByTemperature ddTp_analytical;
  Medium.DerDensityByTemperature ddTp_numerical;
  Medium.DerDensityByPressure ddpT_analytical;
  Medium.DerDensityByPressure ddpT_numerical;

// Enthalpy derivatives
  Medium.SpecificHeatCapacity dhTp_analytical;
  Medium.SpecificHeatCapacity dhTp_numerical;
  Medium.DerEnthalpyByPressure dhpT_analytical;
  Medium.DerEnthalpyByPressure dhpT_numerical;

protected
  Real eps= 1e-5;
  Medium.ThermodynamicState    state=Medium.setState_pTX(p=p, T=T);
  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsThird(T=T, d=state.d, phase=state.phase);

  Medium.ThermodynamicState    p_plus=Medium.setState_pTX(p=p+eps*p, T=T);
  Medium.EoS.HelmholtzDerivs f_p_plus=Medium.EoS.setHelmholtzDerivsThird(T=T, d=p_plus.d, phase=state.phase);
  Medium.ThermodynamicState    p_minus=Medium.setState_pTX(p=p-eps*p, T=T);
  Medium.EoS.HelmholtzDerivs f_p_minus=Medium.EoS.setHelmholtzDerivsThird(T=T, d=p_minus.d, phase=state.phase);

  Medium.ThermodynamicState    T_plus=Medium.setState_pTX(p=p, T=T+eps*T);
  Medium.EoS.HelmholtzDerivs f_T_plus=Medium.EoS.setHelmholtzDerivsThird(T=T_plus.T, d=T_plus.d, phase=state.phase);
  Medium.ThermodynamicState    T_minus=Medium.setState_pTX(p=p, T=T-eps*T);
  Medium.EoS.HelmholtzDerivs f_T_minus=Medium.EoS.setHelmholtzDerivsThird(T=T_minus.T, d=T_minus.d, phase=state.phase);

equation
  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print(" ");

  Modelica.Utilities.Streams.print("Density");
  // check (dd/dT)@p=const
  ddTp_analytical = Medium.density_derT_p(state=state);
  ddTp_numerical = (T_plus.d-T_minus.d)/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const analytical= " + String(ddTp_analytical));
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const  numerical= " + String(ddTp_numerical));
  // check (dd/dp)@T=const
  ddpT_analytical = Medium.density_derp_T(state=state);
  ddpT_numerical = (p_plus.d-p_minus.d)/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const analytical= " + String(ddpT_analytical));
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const  numerical= " + String(ddpT_numerical));

Modelica.Utilities.Streams.print("Enthalpy");
  // check (dh/dT)@p=const
  dhTp_analytical = Medium.specificHeatCapacityCp(state=state);
  dhTp_numerical = (T_plus.h-T_minus.h)/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("  (dh/dT)@p=const analytical= " + String(dhTp_analytical));
  Modelica.Utilities.Streams.print("  (dh/dT)@p=const  numerical= " + String(dhTp_numerical));
  // check (dh/dp)@T=const
  dhpT_analytical = Medium.isothermalThrottlingCoefficient(state=state);
  dhpT_numerical = (p_plus.h-p_minus.h)/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (dh/dp)@T=const analytical= " + String(dhpT_analytical));
  Modelica.Utilities.Streams.print("  (dh/dp)@T=const  numerical= " + String(dhpT_numerical));

annotation (experiment(NumberOfIntervals=1));
end Derivatives_SinglePhase_pT;
