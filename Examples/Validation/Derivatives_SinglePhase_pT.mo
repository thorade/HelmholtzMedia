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

  Medium.Types.Der2DensityByTemperature2 d2dT2p_analytical;
  Medium.Types.Der2DensityByTemperature2 d2dT2p_numerical;
  Medium.Types.Der2DensityByPressure2 d2dp2T_analytical;
  Medium.Types.Der2DensityByPressure2 d2dp2T_numerical;
  Medium.Types.Der2DensityByTemperaturePressure d2dTp_analytical;
  Medium.Types.Der2DensityByTemperaturePressure d2dTp_numerical1;
  Medium.Types.Der2DensityByTemperaturePressure d2dTp_numerical2;

// Enthalpy derivatives
  Medium.Types.DerEnthalpyByTemperature dhTp_analytical;
  Medium.Types.DerEnthalpyByTemperature dhTp_numerical;
  Medium.DerEnthalpyByPressure dhpT_analytical;
  Medium.DerEnthalpyByPressure dhpT_numerical;

protected
  Real eps= 1e-5;
  Medium.ThermodynamicState state=Medium.setState_pTX(p=p, T=T);
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
  // check (d2d/dT2)@p=const
  d2dT2p_analytical = -(Medium.EoS.d2pT2d(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pTd(f))/(Medium.EoS.dpdT(f)^2)
                      +(Medium.EoS.d2pTd(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f))/(Medium.EoS.dpdT(f)^2) * Medium.EoS.dpTd(f)/Medium.EoS.dpdT(f);
  d2dT2p_numerical = (Medium.density_derT_p(T_plus)-Medium.density_derT_p(T_minus))/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const analytical= " + String(d2dT2p_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const  numerical= " + String(d2dT2p_numerical));
  // check (d2d/dp2)@T=const
  d2dp2T_analytical = -Medium.EoS.d2pd2T(f)/Medium.EoS.dpdT(f)^3;
  d2dp2T_numerical = (Medium.density_derp_T(p_plus)-Medium.density_derp_T(p_minus))/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const analytical= " + String(d2dp2T_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const  numerical= " + String(d2dp2T_numerical));
  // check (d2d/dT dp)
  d2dTp_analytical = -(Medium.EoS.d2pTd(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f))/(Medium.EoS.dpdT(f)^3);
  d2dTp_numerical1 = (Medium.density_derp_T(T_plus)-Medium.density_derp_T(T_minus))/(T_plus.T-T_minus.T);
  d2dTp_numerical2 = (Medium.density_derT_p(p_plus)-Medium.density_derT_p(p_minus))/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (d2d/dT dp) analytical= " + String(d2dTp_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dp dT) numerical1= " + String(d2dTp_numerical1));
  Modelica.Utilities.Streams.print("  (d2d/dT dp) numerical2= " + String(d2dTp_numerical2));

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
