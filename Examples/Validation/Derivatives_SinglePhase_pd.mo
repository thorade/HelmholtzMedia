within HelmholtzMedia.Examples.Validation;
model Derivatives_SinglePhase_pd
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Butane;
  // choose p and d that result in single-phase
  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Density d=2;

// Temperature derivatives
  Medium.Types.DerTemperatureByDensity dTdp_analytical;
  Medium.Types.DerTemperatureByDensity dTdp_numerical;
  Medium.Types.DerTemperatureByPressure dTpd_analytical;
  Medium.Types.DerTemperatureByPressure dTpd_numerical;

  Medium.Types.Der2TemperatureByDensity2 d2Td2p_analytical;
  Medium.Types.Der2TemperatureByDensity2 d2Td2p_numerical;
  Medium.Types.Der2TemperatureByPressure2 d2Tp2d_analytical;
  Medium.Types.Der2TemperatureByPressure2 d2Tp2d_numerical;
  Medium.Types.Der2TemperatureByPressureDensity d2Tpd_analytical;
  Medium.Types.Der2TemperatureByPressureDensity d2Tpd_numerical1;
  Medium.Types.Der2TemperatureByPressureDensity d2Tpd_numerical2;

protected
  Real eps= 1e-5;
  Medium.ThermodynamicState state=Medium.setState_pd(p=p, d=d);
  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsThird(T=state.T, d=state.d, phase=state.phase);

  Medium.ThermodynamicState    d_plus=Medium.setState_pd(p=p, d=d+eps*d);
  Medium.EoS.HelmholtzDerivs f_d_plus=Medium.EoS.setHelmholtzDerivsThird(T=d_plus.T, d=d_plus.d, phase=state.phase);
  Medium.ThermodynamicState    d_minus=Medium.setState_pd(p=p, d=d-eps*d);
  Medium.EoS.HelmholtzDerivs f_d_minus=Medium.EoS.setHelmholtzDerivsThird(T=d_minus.T, d=d_minus.d, phase=state.phase);

  Medium.ThermodynamicState    p_plus=Medium.setState_pd(p=p+eps*p, d=d);
  Medium.EoS.HelmholtzDerivs f_p_plus=Medium.EoS.setHelmholtzDerivsThird(T=p_plus.T, d=p_plus.d, phase=state.phase);
  Medium.ThermodynamicState    p_minus=Medium.setState_pd(p=p-eps*p, d=d);
  Medium.EoS.HelmholtzDerivs f_p_minus=Medium.EoS.setHelmholtzDerivsThird(T=p_minus.T, d=p_minus.d, phase=state.phase);

equation
  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print(" ");

  Modelica.Utilities.Streams.print("Temperature");
  // check (dT/dd)@p=const
  dTdp_analytical = -Medium.EoS.dpdT(f)/Medium.EoS.dpTd(f);
  dTdp_numerical = (d_plus.T-d_minus.T)/(d_plus.d-d_minus.d);
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const analytical= " + String(dTdp_analytical));
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const  numerical= " + String(dTdp_numerical));
  // check (dT/dp)@d=const
  dTpd_analytical = 1/Medium.EoS.dpTd(f);
  dTpd_numerical = (p_plus.T-p_minus.T)/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const analytical= " + String(dTpd_analytical));
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const  numerical= " + String(dTpd_numerical));
  // check (d2T/dd2)@p=const
  d2Td2p_analytical = -(Medium.EoS.d2pd2T(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pTd(f))/(Medium.EoS.dpTd(f)^2)
                      +(Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pT2d(f))/(Medium.EoS.dpTd(f)^2) * Medium.EoS.dpdT(f)/Medium.EoS.dpTd(f);
  d2Td2p_numerical = (-Medium.EoS.dpdT(f_d_plus)/Medium.EoS.dpTd(f_d_plus)+Medium.EoS.dpdT(f_d_minus)/Medium.EoS.dpTd(f_d_minus))/(d_plus.d-d_minus.d);
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const analytical= " + String(d2Td2p_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const  numerical= " + String(d2Td2p_numerical));
  // check (d2T/dp2)@d=const
  d2Tp2d_analytical = -Medium.EoS.d2pT2d(f)/Medium.EoS.dpTd(f)^3;
  d2Tp2d_numerical = (1/Medium.EoS.dpTd(f_p_plus)-1/Medium.EoS.dpTd(f_p_minus))/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const analytical= " + String(d2Tp2d_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const  numerical= " + String(d2Tp2d_numerical));
  // check (d2T/dp dd)
  d2Tpd_analytical = -(Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pT2d(f))/(Medium.EoS.dpTd(f)^3);
  d2Tpd_numerical1 = (1/Medium.EoS.dpTd(f_d_plus)-1/Medium.EoS.dpTd(f_d_minus))/(d_plus.d-d_minus.d);
  d2Tpd_numerical2 = (-Medium.EoS.dpdT(f_p_plus)/Medium.EoS.dpTd(f_p_plus)+Medium.EoS.dpdT(f_p_minus)/Medium.EoS.dpTd(f_p_minus))/(p_plus.p-p_minus.p);
  Modelica.Utilities.Streams.print("  (d2T/dp dd) analytical= " + String(d2Tpd_analytical));
  Modelica.Utilities.Streams.print("  (d2T/dp dd) numerical1= " + String(d2Tpd_numerical1));
  Modelica.Utilities.Streams.print("  (d2T/dd dp) numerical2= " + String(d2Tpd_numerical2));

annotation (experiment(NumberOfIntervals=1));
end Derivatives_SinglePhase_pd;
