within HelmholtzMedia.Examples.Validation;
model Derivatives_SinglePhase
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Butane;
  // p and T always result in single-phase
  parameter Medium.AbsolutePressure p=101325;
  parameter Medium.Temperature T=298.15;
  Medium.ThermodynamicState state=Medium.setState_pTX(p=p, T=T);
  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsThird(T=state.T, d=state.d, phase=state.phase);

// Pressure wrt. dT
  Medium.Types.DerPressureByDensity dpdT_analytical;
  Medium.Types.DerPressureByDensity dpdT_numerical;
  Medium.Types.DerPressureByTemperature dpTd_analytical;
  Medium.Types.DerPressureByTemperature dpTd_numerical;
// Pressure wrt. dT 2nd order
  Medium.Types.Der2PressureByDensity2 d2pd2T_analytical;
  Medium.Types.Der2PressureByDensity2 d2pd2T_numerical;
  Medium.Types.Der2PressureByTemperature2 d2pT2d_analytical;
  Medium.Types.Der2PressureByTemperature2 d2pT2d_numerical;
  Medium.Types.Der2PressureByTemperatureDensity d2pTd_analytical;
  Medium.Types.Der2PressureByTemperatureDensity d2pTd_numerical;
// Entropy wrt. dT
  Medium.Types.DerEntropyByDensity dsdT_analytical;
  Medium.Types.DerEntropyByDensity dsdT_numerical;
  Medium.Types.DerEntropyByTemperature dsTd_analytical;
  Medium.Types.DerEntropyByTemperature dsTd_numerical;
// Entropy wrt. dT 2nd order
  Medium.Types.Der2EntropyByDensity2 d2sd2T_analytical;
  Medium.Types.Der2EntropyByDensity2 d2sd2T_numerical;
  Medium.Types.Der2EntropyByTemperature2 d2sT2d_analytical;
  Medium.Types.Der2EntropyByTemperature2 d2sT2d_numerical;
  Medium.Types.Der2EntropyByTemperatureDensity d2sTd_analytical;
  Medium.Types.Der2EntropyByTemperatureDensity d2sTd_numerical;
// Energy wrt. dT
  Medium.Types.DerEnergyByDensity dudT_analytical;
  Medium.Types.DerEnergyByDensity dudT_numerical;
  Medium.Types.DerEnergyByTemperature duTd_analytical;
  Medium.Types.DerEnergyByTemperature duTd_numerical;
// Energy wrt. dT 2nd order
  Medium.Types.Der2EnergyByDensity2 d2ud2T_analytical;
  Medium.Types.Der2EnergyByDensity2 d2ud2T_numerical;
  Medium.Types.Der2EnergyByTemperature2 d2uT2d_analytical;
  Medium.Types.Der2EnergyByTemperature2 d2uT2d_numerical;
  Medium.Types.Der2EnergyByTemperatureDensity d2uTd_analytical;
  Medium.Types.Der2EnergyByTemperatureDensity d2uTd_numerical;
// Enthalpy wrt. dT
  Medium.Types.DerEnthalpyByDensity dhdT_analytical;
  Medium.Types.DerEnthalpyByDensity dhdT_numerical;
  Medium.Types.DerEnthalpyByTemperature dhTd_analytical;
  Medium.Types.DerEnthalpyByTemperature dhTd_numerical;
// Enthalpy wrt. dT 2nd order
  Medium.Types.Der2EnthalpyByDensity2 d2hd2T_analytical;
  Medium.Types.Der2EnthalpyByDensity2 d2hd2T_numerical;
  Medium.Types.Der2EnthalpyByTemperature2 d2hT2d_analytical;
  Medium.Types.Der2EnthalpyByTemperature2 d2hT2d_numerical;
  Medium.Types.Der2EnthalpyByTemperatureDensity d2hTd_analytical;
  Medium.Types.Der2EnthalpyByTemperatureDensity d2hTd_numerical;
// Enthalpy wrt. pT
  Medium.Types.DerEnthalpyByTemperature dhTp_analytical;
  Medium.Types.DerEnthalpyByTemperature dhTp_numerical;
  Medium.DerEnthalpyByPressure dhpT_analytical;
  Medium.DerEnthalpyByPressure dhpT_numerical;
// Gibbs wrt. dT
  Medium.Types.DerEnergyByDensity dgdT_analytical;
  Medium.Types.DerEnergyByDensity dgdT_numerical;
  Medium.Types.DerEnergyByTemperature dgTd_analytical;
  Medium.Types.DerEnergyByTemperature dgTd_numerical;
  Medium.Types.Der2EnergyByDensity2 d2gd2T_analytical;
  Medium.Types.Der2EnergyByDensity2 d2gd2T_numerical;
  Medium.Types.Der2EnergyByTemperature2 d2gT2d_analytical;
  Medium.Types.Der2EnergyByTemperature2 d2gT2d_numerical;
  Medium.Types.Der2EnergyByTemperatureDensity d2gTd_analytical;
  Medium.Types.Der2EnergyByTemperatureDensity d2gTd_numerical;
// Density wrt. pT
  Medium.DerDensityByTemperature ddTp_analytical;
  Medium.DerDensityByTemperature ddTp_numerical;
  Medium.DerDensityByPressure ddpT_analytical;
  Medium.DerDensityByPressure ddpT_numerical;
// Density wrt. pT 2nd order
  Medium.Types.Der2DensityByTemperature2 d2dT2p_analytical;
  Medium.Types.Der2DensityByTemperature2 d2dT2p_numerical;
  Medium.Types.Der2DensityByPressure2 d2dp2T_analytical;
  Medium.Types.Der2DensityByPressure2 d2dp2T_numerical;
  Medium.Types.Der2DensityByTemperaturePressure d2dTp_analytical;
  Medium.Types.Der2DensityByTemperaturePressure d2dTp_numerical1;
  Medium.Types.Der2DensityByTemperaturePressure d2dTp_numerical2;
// Temperature wrt. pd
  Medium.Types.DerTemperatureByDensity dTdp_analytical;
  Medium.Types.DerTemperatureByDensity dTdp_numerical;
  Medium.Types.DerTemperatureByPressure dTpd_analytical;
  Medium.Types.DerTemperatureByPressure dTpd_numerical;
// Temperature wrt. pd 2nd order
  Medium.Types.Der2TemperatureByDensity2 d2Td2p_analytical;
  Medium.Types.Der2TemperatureByDensity2 d2Td2p_numerical;
  Medium.Types.Der2TemperatureByPressure2 d2Tp2d_analytical;
  Medium.Types.Der2TemperatureByPressure2 d2Tp2d_numerical;
  Medium.Types.Der2TemperatureByPressureDensity d2Tpd_analytical;
  Medium.Types.Der2TemperatureByPressureDensity d2Tpd_numerical1;
  Medium.Types.Der2TemperatureByPressureDensity d2Tpd_numerical2;
// Furter derivatives, 1st order
// Further derivatives, 2nd order
  Medium.Types.Der2EnthalpyByTemperature2 dcpTd_analytical;
  Medium.Types.Der2EnthalpyByTemperature2 dcpTd_numerical;
  Medium.Types.Der2EnthalpyByTemperatureDensity dcpdT_analytical;
  Medium.Types.Der2EnthalpyByTemperatureDensity dcpdT_numerical;

protected
  Real eps= 1e-5;

  Medium.ThermodynamicState    dplus_Tconst=Medium.setState_dTX(d=state.d+eps*state.d, T=state.T);
  Medium.EoS.HelmholtzDerivs f_dplus_Tconst=Medium.EoS.setHelmholtzDerivsThird(T=state.T, d=dplus_Tconst.d, phase=state.phase);
  Medium.ThermodynamicState    dminus_Tconst=Medium.setState_dTX(d=state.d-eps*state.d, T=state.T);
  Medium.EoS.HelmholtzDerivs f_dminus_Tconst=Medium.EoS.setHelmholtzDerivsThird(T=state.T, d=dminus_Tconst.d, phase=state.phase);

  Medium.ThermodynamicState    Tplus_dconst=Medium.setState_dTX(d=state.d, T=state.T+eps*state.T);
  Medium.EoS.HelmholtzDerivs f_Tplus_dconst=Medium.EoS.setHelmholtzDerivsThird(T=Tplus_dconst.T, d=state.d, phase=state.phase);
  Medium.ThermodynamicState    Tminus_dconst=Medium.setState_dTX(d=state.d, T=state.T-eps*state.T);
  Medium.EoS.HelmholtzDerivs f_Tminus_dconst=Medium.EoS.setHelmholtzDerivsThird(T=Tminus_dconst.T, d=state.d, phase=state.phase);

  Medium.ThermodynamicState    pplus_Tconst=Medium.setState_pTX(p=state.p+eps*state.p, T=state.T);
  Medium.EoS.HelmholtzDerivs f_pplus_Tconst=Medium.EoS.setHelmholtzDerivsThird(T=pplus_Tconst.T, d=pplus_Tconst.d, phase=state.phase);
  Medium.ThermodynamicState    pminus_Tconst=Medium.setState_pTX(p=state.p-eps*state.p, T=state.T);
  Medium.EoS.HelmholtzDerivs f_pminus_Tconst=Medium.EoS.setHelmholtzDerivsThird(T=pminus_Tconst.T, d=pminus_Tconst.d, phase=state.phase);

  Medium.ThermodynamicState    Tplus_pconst=Medium.setState_pTX(p=state.p, T=state.T+eps*state.T);
  Medium.EoS.HelmholtzDerivs f_Tplus_pconst=Medium.EoS.setHelmholtzDerivsThird(T=Tplus_pconst.T, d=Tplus_pconst.d, phase=state.phase);
  Medium.ThermodynamicState    Tminus_pconst=Medium.setState_pTX(p=state.p, T=state.T-eps*state.T);
  Medium.EoS.HelmholtzDerivs f_Tminus_pconst=Medium.EoS.setHelmholtzDerivsThird(T=Tminus_pconst.T, d=Tminus_pconst.d, phase=state.phase);

  Medium.ThermodynamicState    dplus_pconst=Medium.setState_pd(p=state.p, d=state.d+eps*state.d);
  Medium.EoS.HelmholtzDerivs f_dplus_pconst=Medium.EoS.setHelmholtzDerivsThird(T=dplus_pconst.T, d=dplus_pconst.d, phase=state.phase);
  Medium.ThermodynamicState    dminus_pconst=Medium.setState_pd(p=state.p, d=state.d-eps*state.d);
  Medium.EoS.HelmholtzDerivs f_dminus_pconst=Medium.EoS.setHelmholtzDerivsThird(T=dminus_pconst.T, d=dminus_pconst.d, phase=state.phase);

  Medium.ThermodynamicState    pplus_dconst=Medium.setState_pd(p=state.p+eps*state.p, d=state.d);
  Medium.EoS.HelmholtzDerivs f_pplus_dconst=Medium.EoS.setHelmholtzDerivsThird(T=pplus_dconst.T, d=pplus_dconst.d, phase=state.phase);
  Medium.ThermodynamicState    pminus_dconst=Medium.setState_pd(p=state.p-eps*state.p, d=state.d);
  Medium.EoS.HelmholtzDerivs f_pminus_dconst=Medium.EoS.setHelmholtzDerivsThird(T=pminus_dconst.T, d=pminus_dconst.d, phase=state.phase);

equation
  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Pressure");
  // check (dp/dd)@T=const
  dpdT_analytical = Medium.EoS.dpdT(f);
  dpdT_numerical = (dplus_Tconst.p-dminus_Tconst.p)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (dp/dd)@T=const analytical= " + String(dpdT_analytical));
  Modelica.Utilities.Streams.print("  (dp/dd)@T=const  numerical= " + String(dpdT_numerical));
  // check (dp/dT)@d=const
  dpTd_analytical = Medium.EoS.dpTd(f);
  dpTd_numerical = (Tplus_dconst.p-Tminus_dconst.p)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (dp/dT)@d=const analytical= " + String(dpTd_analytical));
  Modelica.Utilities.Streams.print("  (dp/dT)@d=const  numerical= " + String(dpTd_numerical));
  // check (d2p/dd2)@T=const
  d2pd2T_analytical = Medium.EoS.d2pd2T(f);
  d2pd2T_numerical = (Medium.EoS.dpdT(f_dplus_Tconst)-Medium.EoS.dpdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2p/dd2)@T=const analytical= " + String(d2pd2T_analytical));
  Modelica.Utilities.Streams.print("  (d2p/dd2)@T=const  numerical= " + String(d2pd2T_numerical));
  // check (d2p/dT2)@d=const
  d2pT2d_analytical = Medium.EoS.d2pT2d(f);
  d2pT2d_numerical = (Medium.EoS.dpTd(f_Tplus_dconst)-Medium.EoS.dpTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2p/dT2)@d=const analytical= " + String(d2pT2d_analytical));
  Modelica.Utilities.Streams.print("  (d2p/dT2)@d=const  numerical= " + String(d2pT2d_numerical));
  // check (d2p/dT dd)
  d2pTd_analytical = Medium.EoS.d2pTd(f);
  d2pTd_numerical = (Medium.EoS.dpTd(f_dplus_Tconst)-Medium.EoS.dpTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2p/dTdd) analytical= " + String(d2pTd_analytical));
  Modelica.Utilities.Streams.print("  (d2p/dTdd)  numerical= " + String(d2pTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Entropy");
  // check (ds/dd)@T=const
  dsdT_analytical = Medium.EoS.dsdT(f);
  dsdT_numerical = (dplus_Tconst.s-dminus_Tconst.s)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (ds/dd)@T=const analytical= " + String(dsdT_analytical));
  Modelica.Utilities.Streams.print("  (ds/dd)@T=const  numerical= " + String(dsdT_numerical));
  // check (ds/dT)@d=const
  dsTd_analytical = Medium.EoS.dsTd(f);
  dsTd_numerical = (Tplus_dconst.s-Tminus_dconst.s)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const analytical= " + String(dsTd_analytical));
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const  numerical= " + String(dsTd_numerical));
  // check (d2u/dd2)@T=const
  d2sd2T_analytical = Medium.EoS.d2sd2T(f);
  d2sd2T_numerical = (Medium.EoS.dsdT(f_dplus_Tconst)-Medium.EoS.dsdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2s/dd2)@T=const analytical= " + String(d2sd2T_analytical));
  Modelica.Utilities.Streams.print("  (d2s/dd2)@T=const  numerical= " + String(d2sd2T_numerical));
  // check (d2s/dT2)@d=const
  d2sT2d_analytical = Medium.EoS.d2sT2d(f);
  d2sT2d_numerical = (Medium.EoS.dsTd(f_Tplus_dconst)-Medium.EoS.dsTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2s/dT2)@d=const analytical= " + String(d2sT2d_analytical));
  Modelica.Utilities.Streams.print("  (d2s/dT2)@d=const  numerical= " + String(d2sT2d_numerical));
  // check (d2s/dT dd)
  d2sTd_analytical = Medium.EoS.d2sTd(f);
  d2sTd_numerical = (Medium.EoS.dsTd(f_dplus_Tconst)-Medium.EoS.dsTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2s/dT dd) analytical= " + String(d2sTd_analytical));
  Modelica.Utilities.Streams.print("  (d2s/dT dd)  numerical= " + String(d2sTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("internal Energy");
  // check (du/dd)@T=const
  dudT_analytical = Medium.EoS.dudT(f);
  dudT_numerical = (dplus_Tconst.u-dminus_Tconst.u)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (du/dd)@T=const analytical= " + String(dudT_analytical));
  Modelica.Utilities.Streams.print("  (du/dd)@T=const  numerical= " + String(dudT_numerical));
  // check (du/dT)@d=const
  duTd_analytical = Medium.EoS.duTd(f);
  duTd_numerical = (Tplus_dconst.u-Tminus_dconst.u)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (du/dT)@d=const analytical= " + String(duTd_analytical));
  Modelica.Utilities.Streams.print("  (du/dT)@d=const  numerical= " + String(duTd_numerical));
  // check (d2u/dd2)@T=const
  d2ud2T_analytical = Medium.EoS.d2ud2T(f);
  d2ud2T_numerical = (Medium.EoS.dudT(f_dplus_Tconst)-Medium.EoS.dudT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2u/dd2)@T=const analytical= " + String(d2ud2T_analytical));
  Modelica.Utilities.Streams.print("  (d2u/dd2)@T=const  numerical= " + String(d2ud2T_numerical));
  // check (d2u/dT2)@d=const
  d2uT2d_analytical = Medium.EoS.d2uT2d(f);
  d2uT2d_numerical = (Medium.EoS.duTd(f_Tplus_dconst)-Medium.EoS.duTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2u/dT2)@d=const analytical= " + String(d2uT2d_analytical));
  Modelica.Utilities.Streams.print("  (d2u/dT2)@d=const  numerical= " + String(d2uT2d_numerical));
  // check (d2u/dT dd)
  d2uTd_analytical = Medium.EoS.d2uTd(f);
  d2uTd_numerical = (Medium.EoS.duTd(f_dplus_Tconst)-Medium.EoS.duTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2u/dT dd) analytical= " + String(d2uTd_analytical));
  Modelica.Utilities.Streams.print("  (d2u/dT dd)  numerical= " + String(d2uTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Enthalpy");
  // check (dh/dd)@T=const
  dhdT_analytical = Medium.EoS.dhdT(f);
  dhdT_numerical = (dplus_Tconst.h-dminus_Tconst.h)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const analytical= " + String(dhdT_analytical));
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const  numerical= " + String(dhdT_numerical));
  // check (dh/dT)@d=const
  dhTd_analytical = Medium.EoS.dhTd(f);
  dhTd_numerical = (Tplus_dconst.h-Tminus_dconst.h)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const analytical= " + String(dhTd_analytical));
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const  numerical= " + String(dhTd_numerical));
  // check (d2h/dd2)@T=const
  d2hd2T_analytical = Medium.EoS.d2hd2T(f);
  d2hd2T_numerical = (Medium.EoS.dhdT(f_dplus_Tconst)-Medium.EoS.dhdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2h/dd2)@T=const analytical= " + String(d2hd2T_analytical));
  Modelica.Utilities.Streams.print("  (d2h/dd2)@T=const  numerical= " + String(d2hd2T_numerical));
  // check (d2h/dT2)@d=const
  d2hT2d_analytical = Medium.EoS.d2hT2d(f);
  d2hT2d_numerical = (Medium.EoS.dhTd(f_Tplus_dconst)-Medium.EoS.dhTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2h/dT2)@d=const analytical= " + String(d2hT2d_analytical));
  Modelica.Utilities.Streams.print("  (d2h/dT2)@d=const  numerical= " + String(d2hT2d_numerical));
  // check (d2h/dT dd)
  d2hTd_analytical = Medium.EoS.d2hTd(f);
  d2hTd_numerical = (Medium.EoS.dhTd(f_dplus_Tconst)-Medium.EoS.dhTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2h/dT dd) analytical= " + String(d2hTd_analytical));
  Modelica.Utilities.Streams.print("  (d2h/dT dd)  numerical= " + String(d2hTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Gibbs energy");
  // check (dg/dd)@T=const
  dgdT_analytical = Medium.EoS.dgdT(f);
  dgdT_numerical = ((dplus_Tconst.h-dplus_Tconst.T*dplus_Tconst.s)-(dminus_Tconst.h-dminus_Tconst.T*dminus_Tconst.s))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (dg/dd)@T=const analytical= " + String(dgdT_analytical));
  Modelica.Utilities.Streams.print("  (dg/dd)@T=const  numerical= " + String(dgdT_numerical));
  // check (dg/dT)@d=const
  dgTd_analytical = Medium.EoS.dgTd(f);
  dgTd_numerical = ((Tplus_dconst.h-Tplus_dconst.T*Tplus_dconst.s)-(Tminus_dconst.h-Tminus_dconst.T*Tminus_dconst.s))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (dg/dT)@d=const analytical= " + String(dgTd_analytical));
  Modelica.Utilities.Streams.print("  (dg/dT)@d=const  numerical= " + String(dgTd_numerical));
  // check (d2g/dd2)@T=const
  d2gd2T_analytical = Medium.EoS.d2gd2T(f);
  d2gd2T_numerical = (Medium.EoS.dgdT(f_dplus_Tconst)-Medium.EoS.dgdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2g/dd2)@T=const analytical= " + String(d2gd2T_analytical));
  Modelica.Utilities.Streams.print("  (d2g/dd2)@T=const  numerical= " + String(d2gd2T_numerical));
  // check (d2g/dT2)@d=const
  d2gT2d_analytical = Medium.EoS.d2gT2d(f);
  d2gT2d_numerical = (Medium.EoS.dgTd(f_Tplus_dconst)-Medium.EoS.dgTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2g/dT2)@d=const analytical= " + String(d2gT2d_analytical));
  Modelica.Utilities.Streams.print("  (d2g/dT2)@d=const  numerical= " + String(d2gT2d_numerical));
  // check (d2g/dT dd)
  d2gTd_analytical = Medium.EoS.d2gTd(f);
  d2gTd_numerical = (Medium.EoS.dgTd(f_dplus_Tconst)-Medium.EoS.dgTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2g/dT dd) analytical= " + String(d2gTd_analytical));
  Modelica.Utilities.Streams.print("  (d2g/dT dd)  numerical= " + String(d2gTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Density");
  // check (dd/dT)@p=const
  ddTp_analytical = Medium.density_derT_p(state=state);
  ddTp_numerical = (Tplus_pconst.d-Tminus_pconst.d)/(Tplus_pconst.T-Tminus_pconst.T);
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const analytical= " + String(ddTp_analytical));
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const  numerical= " + String(ddTp_numerical));
  // check (dd/dp)@T=const
  ddpT_analytical = Medium.density_derp_T(state=state);
  ddpT_numerical = (pplus_Tconst.d-pminus_Tconst.d)/(pplus_Tconst.p-pminus_Tconst.p);
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const analytical= " + String(ddpT_analytical));
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const  numerical= " + String(ddpT_numerical));
  // check (d2d/dT2)@p=const
  d2dT2p_analytical = -(Medium.EoS.d2pT2d(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pTd(f))/(Medium.EoS.dpdT(f)^2)
                      +(Medium.EoS.d2pTd(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f))/(Medium.EoS.dpdT(f)^2) * Medium.EoS.dpTd(f)/Medium.EoS.dpdT(f);
  d2dT2p_numerical = (Medium.density_derT_p(Tplus_pconst)-Medium.density_derT_p(Tminus_pconst))/(Tplus_pconst.T-Tminus_pconst.T);
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const analytical= " + String(d2dT2p_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const  numerical= " + String(d2dT2p_numerical));
  // check (d2d/dp2)@T=const
  d2dp2T_analytical = -Medium.EoS.d2pd2T(f)/Medium.EoS.dpdT(f)^3;
  d2dp2T_numerical = (Medium.density_derp_T(pplus_Tconst)-Medium.density_derp_T(pminus_Tconst))/(pplus_Tconst.p-pminus_Tconst.p);
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const analytical= " + String(d2dp2T_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const  numerical= " + String(d2dp2T_numerical));
  // check (d2d/dT dp)
  d2dTp_analytical = -(Medium.EoS.d2pTd(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f))/(Medium.EoS.dpdT(f)^3);
  d2dTp_numerical1 = (Medium.density_derp_T(Tplus_pconst)-Medium.density_derp_T(Tminus_pconst))/(Tplus_pconst.T-Tminus_pconst.T);
  d2dTp_numerical2 = (Medium.density_derT_p(pplus_Tconst)-Medium.density_derT_p(pminus_Tconst))/(pplus_Tconst.p-pminus_Tconst.p);
  Modelica.Utilities.Streams.print("  (d2d/dT dp) analytical= " + String(d2dTp_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dp dT) numerical1= " + String(d2dTp_numerical1));
  Modelica.Utilities.Streams.print("  (d2d/dT dp) numerical2= " + String(d2dTp_numerical2));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Temperature");
  // check (dT/dd)@p=const
  dTdp_analytical = -Medium.EoS.dpdT(f)/Medium.EoS.dpTd(f);
  dTdp_numerical = (dplus_pconst.T-dminus_pconst.T)/(dplus_pconst.d-dminus_pconst.d);
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const analytical= " + String(dTdp_analytical));
  Modelica.Utilities.Streams.print("  (dd/dT)@p=const  numerical= " + String(dTdp_numerical));
  // check (dT/dp)@d=const
  dTpd_analytical = 1/Medium.EoS.dpTd(f);
  dTpd_numerical = (pplus_dconst.T-pminus_dconst.T)/(pplus_dconst.p-pminus_dconst.p);
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const analytical= " + String(dTpd_analytical));
  Modelica.Utilities.Streams.print("  (dd/dp)@T=const  numerical= " + String(dTpd_numerical));
  // check (d2T/dd2)@p=const
  d2Td2p_analytical = -(Medium.EoS.d2pd2T(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pTd(f))/(Medium.EoS.dpTd(f)^2)
                      +(Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pT2d(f))/(Medium.EoS.dpTd(f)^2) * Medium.EoS.dpdT(f)/Medium.EoS.dpTd(f);
  d2Td2p_numerical = (-Medium.EoS.dpdT(f_dplus_pconst)/Medium.EoS.dpTd(f_dplus_pconst)+Medium.EoS.dpdT(f_dminus_pconst)/Medium.EoS.dpTd(f_dminus_pconst))/(dplus_pconst.d-dminus_pconst.d);
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const analytical= " + String(d2Td2p_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dT2)@p=const  numerical= " + String(d2Td2p_numerical));
  // check (d2T/dp2)@d=const
  d2Tp2d_analytical = -Medium.EoS.d2pT2d(f)/Medium.EoS.dpTd(f)^3;
  d2Tp2d_numerical = (1/Medium.EoS.dpTd(f_pplus_dconst)-1/Medium.EoS.dpTd(f_pminus_dconst))/(pplus_dconst.p-pminus_dconst.p);
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const analytical= " + String(d2Tp2d_analytical));
  Modelica.Utilities.Streams.print("  (d2d/dp2)@T=const  numerical= " + String(d2Tp2d_numerical));
  // check (d2T/dp dd)
  d2Tpd_analytical = -(Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pT2d(f))/(Medium.EoS.dpTd(f)^3);
  d2Tpd_numerical1 = (1/Medium.EoS.dpTd(f_dplus_pconst)-1/Medium.EoS.dpTd(f_dminus_pconst))/(dplus_pconst.d-dminus_pconst.d);
  d2Tpd_numerical2 = (-Medium.EoS.dpdT(f_pplus_dconst)/Medium.EoS.dpTd(f_pplus_dconst)+Medium.EoS.dpdT(f_pminus_dconst)/Medium.EoS.dpTd(f_pminus_dconst))/(pplus_dconst.p-pminus_dconst.p);
  Modelica.Utilities.Streams.print("  (d2T/dp dd) analytical= " + String(d2Tpd_analytical));
  Modelica.Utilities.Streams.print("  (d2T/dp dd) numerical1= " + String(d2Tpd_numerical1));
  Modelica.Utilities.Streams.print("  (d2T/dd dp) numerical2= " + String(d2Tpd_numerical2));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Further first derivatives");
  // check (dh/dT)@p=const
  dhTp_analytical = Medium.specificHeatCapacityCp(state=state);
  dhTp_numerical = (Tplus_pconst.h-Tminus_pconst.h)/(Tplus_pconst.T-Tminus_pconst.T);
  Modelica.Utilities.Streams.print("  (dh/dT)@p=const analytical= " + String(dhTp_analytical));
  Modelica.Utilities.Streams.print("  (dh/dT)@p=const  numerical= " + String(dhTp_numerical));
  // check (dh/dp)@T=const
  dhpT_analytical = Medium.isothermalThrottlingCoefficient(state=state);
  dhpT_numerical = (pplus_Tconst.h-pminus_Tconst.h)/(pplus_Tconst.p-pminus_Tconst.p);
  Modelica.Utilities.Streams.print("  (dh/dp)@T=const analytical= " + String(dhpT_analytical));
  Modelica.Utilities.Streams.print("  (dh/dp)@T=const  numerical= " + String(dhpT_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Further second derivatives");
  // check (d cp/dT)@d=const
  dcpTd_analytical =  Medium.EoS.d2hT2d(f) + Medium.EoS.dhdT(f)*Medium.EoS.dpTd(f)*Medium.EoS.d2pTd(f)/Medium.EoS.dpdT(f)^2
                   - (Medium.EoS.d2hTd(f)*Medium.EoS.dpTd(f) + Medium.EoS.dhdT(f)*Medium.EoS.d2pT2d(f))/Medium.EoS.dpdT(f);
  dcpTd_numerical = (Medium.specificHeatCapacityCp(Tplus_dconst)-Medium.specificHeatCapacityCp(Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d cp/dT)@d=const analytical= " + String(dcpTd_analytical));
  Modelica.Utilities.Streams.print("  (d cp/dT)@d=const  numerical= " + String(dcpTd_numerical));
  // check (d cp/dd)@T=const
  dcpdT_analytical =  Medium.EoS.d2hTd(f) + Medium.EoS.dhdT(f)*Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f)/Medium.EoS.dpdT(f)^2
                   - (Medium.EoS.d2hd2T(f)*Medium.EoS.dpTd(f) + Medium.EoS.dhdT(f)*Medium.EoS.d2pTd(f))/Medium.EoS.dpdT(f);
  dcpdT_numerical = (Medium.specificHeatCapacityCp(dplus_Tconst)-Medium.specificHeatCapacityCp(dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d cp/dd)@T=const analytical= " + String(dcpdT_analytical));
  Modelica.Utilities.Streams.print("  (d cp/dd)@T=const  numerical= " + String(dcpdT_numerical));

annotation (experiment(NumberOfIntervals=1));
end Derivatives_SinglePhase;
