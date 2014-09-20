within HelmholtzMedia.Examples.Validation;
model Derivatives_SinglePhase
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Helium;
  // p and T always result in single-phase
  parameter Medium.Temperature T=400;
  Medium.AbsolutePressure p=101325;
  Medium.ThermodynamicState state=Medium.setState_pTX(p=p, T=T);
  //Medium.Density d=Medium.dewDensity(sat=Medium.setSat_T(T));
  //Medium.ThermodynamicState state=Medium.setState_dTX(d=d, T=T);
  Medium.EoS.HelmholtzDerivs f=Medium.EoS.setHelmholtzDerivsThird(T=state.T, d=state.d, phase=state.phase);

// Pressure wrt. dT
  Medium.DerPressureByDensity dpdT_analytical1;
  Medium.DerPressureByDensity dpdT_analytical2;
  Medium.DerPressureByDensity dpdT_numerical;
  Medium.DerPressureByTemperature dpTd_analytical;
  Medium.DerPressureByTemperature dpTd_numerical;
// Pressure wrt. dT 2nd order
  Medium.Der2PressureByDensity2 d2pd2T_analytical1;
  Medium.Der2PressureByDensity2 d2pd2T_analytical2;
  Medium.Der2PressureByDensity2 d2pd2T_numerical;
  Medium.Der2PressureByDensity2 d2pd2s_analytical1;
  Medium.Der2PressureByDensity2 d2pd2s_analytical2;
  Medium.Der2PressureByDensity2 d2pd2s_analytical3;
  Medium.Der2PressureByVolume2 d2pv2s_analytical1;
  Medium.Der2PressureByVolume2 d2pv2s_analytical2;
  Medium.Der2PressureByVolume2 d2pv2s_analytical3;
  Real fd_gd_analytical1_d;
  Real fd_gd_analytical2_d;
  Real fd_gd_analytical1_v;
  Real fd_gd_analytical2_v;
  Real fd_gd_analytical3_v;
  Medium.Der2PressureByTemperature2 d2pT2d_analytical;
  Medium.Der2PressureByTemperature2 d2pT2d_numerical;
  Medium.Der2PressureByTemperatureDensity d2pTd_analytical1;
  Medium.Der2PressureByTemperatureDensity d2pTd_analytical2;
  Medium.Der2PressureByTemperatureDensity d2pTd_numerical;
// Entropy wrt. dT
  Medium.DerEntropyByDensity dsdT_analytical1;
  Medium.DerEntropyByDensity dsdT_analytical2;
  Medium.DerEntropyByDensity dsdT_numerical;
  Medium.DerEntropyByTemperature dsTd_analytical1;
  Medium.DerEntropyByTemperature dsTd_analytical2;
  Medium.DerEntropyByTemperature dsTd_numerical;
// Entropy wrt. dT 2nd order
  Medium.Der2EntropyByDensity2 d2sd2T_analytical1;
  Medium.Der2EntropyByDensity2 d2sd2T_analytical2;
  Medium.Der2EntropyByDensity2 d2sd2T_numerical;
  Medium.Der2EntropyByTemperature2 d2sT2d_analytical1;
  Medium.Der2EntropyByTemperature2 d2sT2d_analytical2;
  Medium.Der2EntropyByTemperature2 d2sT2d_numerical;
  Medium.Der2EntropyByTemperatureDensity d2sTd_analytical1;
  Medium.Der2EntropyByTemperatureDensity d2sTd_analytical2;
  Medium.Der2EntropyByTemperatureDensity d2sTd_numerical;
// Energy wrt. dT
  Medium.DerEnergyByDensity dudT_analytical1;
  Medium.DerEnergyByDensity dudT_analytical2;
  Medium.DerEnergyByDensity dudT_numerical;
  Medium.DerEnergyByTemperature duTd_analytical;
  Medium.DerEnergyByTemperature duTd_numerical;
// Energy wrt. dT 2nd order
  Medium.Der2EnergyByDensity2 d2ud2T_analytical1;
  Medium.Der2EnergyByDensity2 d2ud2T_analytical2;
  Medium.Der2EnergyByDensity2 d2ud2T_numerical;
  Medium.Der2EnergyByTemperature2 d2uT2d_analytical1;
  Medium.Der2EnergyByTemperature2 d2uT2d_analytical2;
  Medium.Der2EnergyByTemperature2 d2uT2d_numerical;
  Medium.Der2EnergyByTemperatureDensity d2uTd_analytical1;
  Medium.Der2EnergyByTemperatureDensity d2uTd_analytical2;
  Medium.Der2EnergyByTemperatureDensity d2uTd_numerical;
// Enthalpy wrt. dT
  Medium.DerEnthalpyByDensity dhdT_analytical1;
  Medium.DerEnthalpyByDensity dhdT_analytical2;
  Medium.DerEnthalpyByDensity dhdT_numerical;
  Medium.DerEnthalpyByTemperature dhTd_analytical1;
  Medium.DerEnthalpyByTemperature dhTd_analytical2;
  Medium.DerEnthalpyByTemperature dhTd_numerical;
// Enthalpy wrt. dT 2nd order
  Medium.Der2EnthalpyByDensity2 d2hd2T_analytical1;
  Medium.Der2EnthalpyByDensity2 d2hd2T_analytical2;
  Medium.Der2EnthalpyByDensity2 d2hd2T_numerical;
  Medium.Der2EnthalpyByTemperature2 d2hT2d_analytical1;
  Medium.Der2EnthalpyByTemperature2 d2hT2d_analytical2;
  Medium.Der2EnthalpyByTemperature2 d2hT2d_numerical;
  Medium.Der2EnthalpyByTemperatureDensity d2hTd_analytical1;
  Medium.Der2EnthalpyByTemperatureDensity d2hTd_analytical2;
  Medium.Der2EnthalpyByTemperatureDensity d2hTd_numerical;
// Enthalpy wrt. pT
  Medium.DerEnthalpyByTemperature dhTp_analytical;
  Medium.DerEnthalpyByTemperature dhTp_numerical;
  Medium.DerEnthalpyByPressure dhpT_analytical;
  Medium.DerEnthalpyByPressure dhpT_numerical;
// Gibbs wrt. dT
  Medium.DerEnergyByDensity dgdT_analytical1;
  Medium.DerEnergyByDensity dgdT_analytical2;
  Medium.DerEnergyByDensity dgdT_numerical;
  Medium.DerEnergyByTemperature dgTd_analytical1;
  Medium.DerEnergyByTemperature dgTd_analytical2;
  Medium.DerEnergyByTemperature dgTd_numerical;
  Medium.Der2EnergyByDensity2 d2gd2T_analytical1;
  Medium.Der2EnergyByDensity2 d2gd2T_analytical2;
  Medium.Der2EnergyByDensity2 d2gd2T_numerical;
  Medium.Der2EnergyByTemperature2 d2gT2d_analytical1;
  Medium.Der2EnergyByTemperature2 d2gT2d_analytical2;
  Medium.Der2EnergyByTemperature2 d2gT2d_numerical;
  Medium.Der2EnergyByTemperatureDensity d2gTd_analytical1;
  Medium.Der2EnergyByTemperatureDensity d2gTd_analytical2;
  Medium.Der2EnergyByTemperatureDensity d2gTd_numerical;
// Density wrt. pT
  Medium.DerDensityByTemperature ddTp_analytical;
  Medium.DerDensityByTemperature ddTp_numerical;
  Medium.DerDensityByPressure ddpT_analytical;
  Medium.DerDensityByPressure ddpT_numerical;
// Density wrt. ph
  Medium.DerDensityByPressure ddph_analytical;
  Medium.DerDensityByPressure ddph_numerical;
  Medium.DerDensityByEnthalpy ddhp_analytical;
  Medium.DerDensityByEnthalpy ddhp_numerical;
// Density wrt. pT 2nd order
  Medium.Der2DensityByTemperature2 d2dT2p_analytical;
  Medium.Der2DensityByTemperature2 d2dT2p_numerical;
  Medium.Der2DensityByPressure2 d2dp2T_analytical;
  Medium.Der2DensityByPressure2 d2dp2T_numerical;
  Medium.Der2DensityByTemperaturePressure d2dTp_analytical;
  Medium.Der2DensityByTemperaturePressure d2dTp_numerical1;
  Medium.Der2DensityByTemperaturePressure d2dTp_numerical2;
// Temperature wrt. pd
  Medium.DerTemperatureByDensity dTdp_analytical;
  Medium.DerTemperatureByDensity dTdp_numerical;
  Medium.DerTemperatureByPressure dTpd_analytical;
  Medium.DerTemperatureByPressure dTpd_numerical;
// Temperature wrt. ph
  Medium.DerTemperatureByPressure dTph_analytical;
  Medium.DerTemperatureByPressure dTph_numerical;
// Temperature wrt. pd 2nd order
  Medium.Der2TemperatureByDensity2 d2Td2p_analytical;
  Medium.Der2TemperatureByDensity2 d2Td2p_numerical;
  Medium.Der2TemperatureByPressure2 d2Tp2d_analytical;
  Medium.Der2TemperatureByPressure2 d2Tp2d_numerical;
  Medium.Der2TemperatureByPressureDensity d2Tpd_analytical;
  Medium.Der2TemperatureByPressureDensity d2Tpd_numerical1;
  Medium.Der2TemperatureByPressureDensity d2Tpd_numerical2;
// Second order derivatives wrt. mixed properties
  Medium.Der2EnthalpyByTemperature2 dcpTd_analytical;
  //Medium.Types.Der2EnthalpyByTemperature2 dcpTd_analytical2;
  Medium.Der2EnthalpyByTemperature2 dcpTd_numerical;
  Medium.Der2EnthalpyByTemperatureDensity dcpdT_analytical;
  Medium.Der2EnthalpyByTemperatureDensity dcpdT_numerical;
  Medium.Der2EnthalpyByTemperaturePressure dcpph_analytical;
  Medium.Der2EnthalpyByTemperaturePressure dcpph_numerical;

  Medium.DerFractionByPressure dmuTd_analytical1;
  Medium.DerFractionByPressure dmuTd_analytical2;
  Medium.DerFractionByPressure dmuTd_numerical;
  Medium.Der2TemperatureByPressureDensity dmudT_analytical1;
  Medium.Der2TemperatureByPressureDensity dmudT_analytical2;
  Medium.Der2TemperatureByPressureDensity dmudT_numerical;
  Medium.Der2TemperatureByPressure2 dmuph_analytical1;
  Medium.Der2TemperatureByPressure2 dmuph_analytical2;
  Medium.Der2TemperatureByPressure2 dmuph_numerical;

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

  Medium.ThermodynamicState    hplus_pconst=Medium.setState_ph(p=state.p, h=state.h+abs(eps*state.h));
  Medium.EoS.HelmholtzDerivs f_hplus_pconst=Medium.EoS.setHelmholtzDerivsThird(T=dplus_pconst.T, d=dplus_pconst.d, phase=state.phase);
  Medium.ThermodynamicState    hminus_pconst=Medium.setState_ph(p=state.p, h=state.h-abs(eps*state.h));
  Medium.EoS.HelmholtzDerivs f_hminus_pconst=Medium.EoS.setHelmholtzDerivsThird(T=dminus_pconst.T, d=dminus_pconst.d, phase=state.phase);

  Medium.ThermodynamicState    pplus_hconst=Medium.setState_ph(p=state.p+eps*state.p, h=state.h);
  Medium.EoS.HelmholtzDerivs f_pplus_hconst=Medium.EoS.setHelmholtzDerivsThird(T=pplus_dconst.T, d=pplus_dconst.d, phase=state.phase);
  Medium.ThermodynamicState    pminus_hconst=Medium.setState_ph(p=state.p-eps*state.p, h=state.h);
  Medium.EoS.HelmholtzDerivs f_pminus_hconst=Medium.EoS.setHelmholtzDerivsThird(T=pminus_dconst.T, d=pminus_dconst.d, phase=state.phase);

equation
  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Pressure wrt. dT");
  // check (dp/dd)@T=const
  dpdT_analytical1 = Medium.EoS.dpdT(f);
  dpdT_analytical2 = -1/state.d^2*Medium.EoS.dpvT(f);
  dpdT_numerical = (dplus_Tconst.p-dminus_Tconst.p)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (dp/dd)@T=const analytical1= " + String(dpdT_analytical1));
  Modelica.Utilities.Streams.print("  (dp/dd)@T=const analytical2= " + String(dpdT_analytical2));
  Modelica.Utilities.Streams.print("  (dp/dd)@T=const   numerical= " + String(dpdT_numerical));
  // check (dp/dT)@d=const
  dpTd_analytical = Medium.EoS.dpTd(f);
  dpTd_numerical = (Tplus_dconst.p-Tminus_dconst.p)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (dp/dT)@d=const analytical= " + String(dpTd_analytical));
  Modelica.Utilities.Streams.print("  (dp/dT)@d=const  numerical= " + String(dpTd_numerical));
  // check (d2p/dd2)@T=const
  d2pd2T_analytical1 = Medium.EoS.d2pd2T(f);
  d2pd2T_analytical2 = Medium.EoS.dpvT(f)*2/state.d^3 + Medium.EoS.d2pv2T(f)/state.d^4;
  d2pd2T_numerical = (Medium.EoS.dpdT(f_dplus_Tconst)-Medium.EoS.dpdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2p/dd2)@T=const analytical1= " + String(d2pd2T_analytical1));
  Modelica.Utilities.Streams.print("  (d2p/dd2)@T=const analytical2= " + String(d2pd2T_analytical2));
  Modelica.Utilities.Streams.print("  (d2p/dd2)@T=const   numerical= " + String(d2pd2T_numerical));
  // check (d2p/dT2)@d=const
  d2pT2d_analytical = Medium.EoS.d2pT2d(f);
  d2pT2d_numerical = (Medium.EoS.dpTd(f_Tplus_dconst)-Medium.EoS.dpTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2p/dT2)@d=const analytical= " + String(d2pT2d_analytical));
  Modelica.Utilities.Streams.print("  (d2p/dT2)@d=const  numerical= " + String(d2pT2d_numerical));
  // check (d2p/dT dd)
  d2pTd_analytical1 = Medium.EoS.d2pTd(f);
  d2pTd_analytical2 = -1/state.d^2*Medium.EoS.d2pTv(f);
  d2pTd_numerical = (Medium.EoS.dpTd(f_dplus_Tconst)-Medium.EoS.dpTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2p/dT dd) analytical1= " + String(d2pTd_analytical1));
  Modelica.Utilities.Streams.print("  (d2p/dT dd) analytical2= " + String(d2pTd_analytical2));
  Modelica.Utilities.Streams.print("  (d2p/dT dd)   numerical= " + String(d2pTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  // calculate (d2p/dd2)@s=const  and  (d2p/dv2)@s=const
  d2pd2s_analytical1 = Medium.EoS.d2pd2T(f)
    - (Medium.EoS.dpTd(f)*Medium.EoS.d2sd2T(f)   + 2*Medium.EoS.d2pTd(f)*Medium.EoS.dsdT(f))/Medium.EoS.dsTd(f)
    + (Medium.EoS.d2pT2d(f)*Medium.EoS.dsdT(f)^2 + 2*Medium.EoS.dpTd(f)*Medium.EoS.dsdT(f)*Medium.EoS.d2sTd(f))/Medium.EoS.dsTd(f)^2
    - (Medium.EoS.dpTd(f)*Medium.EoS.dsdT(f)^2*Medium.EoS.d2sT2d(f))/Medium.EoS.dsTd(f)^3;
  d2pd2s_analytical2 = Medium.EoS.d2pd2T(f)
    - state.T  /state.d^2 *(2/state.d*Medium.EoS.dpTd(f)^2 - 3*Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f))/Medium.EoS.duTd(f)
    + state.T  /state.d^4 *(3*state.T*Medium.EoS.d2pT2d(f)*Medium.EoS.dpTd(f)^2 + Medium.EoS.dpTd(f)^3)/Medium.EoS.duTd(f)^2
    - state.T^2/state.d^4 *(Medium.EoS.d2uT2d(f)*Medium.EoS.dpTd(f)^3)/Medium.EoS.duTd(f)^3;
  d2pv2s_analytical1 = Medium.EoS.d2pv2T(f)
    -3/Medium.EoS.dsTd(f)*Medium.EoS.dpTd(f)*Medium.EoS.d2pTv(f)
    +(Medium.EoS.dpTd(f)/Medium.EoS.dsTd(f))^2*(3*Medium.EoS.d2pT2d(f)+1/state.T*Medium.EoS.dpTd(f)*(1-Medium.EoS.d2uT2d(f)/Medium.EoS.dsTd(f)));
  d2pv2s_analytical2 = Medium.EoS.d2pv2T(f)
    - 3*Medium.EoS.d2pTv(f)*Medium.EoS.dpTd(f)/Medium.EoS.dsTd(f)
    + 3*Medium.EoS.dpTd(f)^2*Medium.EoS.d2pT2d(f)/Medium.EoS.dsTd(f)^2
    -   Medium.EoS.dpTd(f)^3*Medium.EoS.d2sT2d(f)/Medium.EoS.dsTd(f)^3;
  // calculate using second way
  d2pd2s_analytical3 = 2/state.d^3*(-state.d^2)*Medium.velocityOfSound(state)^2 + 1/state.d^4*d2pv2s_analytical1;
  d2pv2s_analytical3 = 2*state.d^3*Medium.velocityOfSound(state)^2 +   state.d^4*d2pd2s_analytical1;
  Modelica.Utilities.Streams.print("  d2pd2s_analytical1 = " + String(d2pd2s_analytical1));
  Modelica.Utilities.Streams.print("  d2pd2s_analytical2 = " + String(d2pd2s_analytical2));
  Modelica.Utilities.Streams.print("  d2pd2s_analytical3 = " + String(d2pd2s_analytical3));
  Modelica.Utilities.Streams.print("  d2pv2s_analytical1 = " + String(d2pv2s_analytical1));
  Modelica.Utilities.Streams.print("  d2pv2s_analytical2 = " + String(d2pv2s_analytical2));
  Modelica.Utilities.Streams.print("  d2pv2s_analytical3 = " + String(d2pv2s_analytical3));
  // check fundamental derivative of gasdynamics
  fd_gd_analytical1_d = 1 + state.d/2/Medium.velocityOfSound(state)^2*d2pd2s_analytical1;
  fd_gd_analytical2_d = 1 + state.d/2/Medium.velocityOfSound(state)^2*d2pd2s_analytical2;
  fd_gd_analytical1_v = 1/state.d^3/2/Medium.velocityOfSound(state)^2*d2pv2s_analytical1;
  fd_gd_analytical2_v = 1/state.d^3/2/Medium.velocityOfSound(state)^2*d2pv2s_analytical2;
  fd_gd_analytical3_v = 1/state.d^3/2/Medium.velocityOfSound(state)^2*d2pv2s_analytical3;
  Modelica.Utilities.Streams.print("  fd_gd_analytical1_d = " + String(fd_gd_analytical1_d));
  Modelica.Utilities.Streams.print("  fd_gd_analytical2_d = " + String(fd_gd_analytical2_d));
  Modelica.Utilities.Streams.print("  fd_gd_analytical1_v = " + String(fd_gd_analytical1_v));
  Modelica.Utilities.Streams.print("  fd_gd_analytical2_v = " + String(fd_gd_analytical2_v));
  //Modelica.Utilities.Streams.print("  fd_gd_analytical3_v = " + String(fd_gd_analytical3_v));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Entropy wrt. dT");
  // check (ds/dd)@T=const
  dsdT_analytical1 = Medium.EoS.dsdT(f);
  dsdT_analytical2 = -1/state.d^2*Medium.EoS.dpTd(f);
  dsdT_numerical = (dplus_Tconst.s-dminus_Tconst.s)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (ds/dd)@T=const analytical1= " + String(dsdT_analytical1));
  Modelica.Utilities.Streams.print("  (ds/dd)@T=const analytical1= " + String(dsdT_analytical2));
  Modelica.Utilities.Streams.print("  (ds/dd)@T=const   numerical= " + String(dsdT_numerical));
  // check (ds/dT)@d=const
  dsTd_analytical1 = Medium.EoS.dsTd(f);
  dsTd_analytical2 = 1/T*Medium.EoS.duTd(f);
  dsTd_numerical = (Tplus_dconst.s-Tminus_dconst.s)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const analytical1= " + String(dsTd_analytical1));
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const analytical2= " + String(dsTd_analytical2));
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const   numerical= " + String(dsTd_numerical));
  // check (d2s/dd2)@T=const
  d2sd2T_analytical1 = Medium.EoS.d2sd2T(f);
  d2sd2T_analytical2 = -1/state.d^2*Medium.EoS.d2pTd(f) + 2/state.d^3*Medium.EoS.dpTd(f);
  d2sd2T_numerical = (Medium.EoS.dsdT(f_dplus_Tconst)-Medium.EoS.dsdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2s/dd2)@T=const analytical1= " + String(d2sd2T_analytical1));
  Modelica.Utilities.Streams.print("  (d2s/dd2)@T=const analytical2= " + String(d2sd2T_analytical2));
  Modelica.Utilities.Streams.print("  (d2s/dd2)@T=const   numerical= " + String(d2sd2T_numerical));
  // check (d2s/dT2)@d=const
  d2sT2d_analytical1 = Medium.EoS.d2sT2d(f);
  d2sT2d_analytical2 = 1/state.T*Medium.EoS.d2uT2d(f) - 1/state.T^2*Medium.EoS.duTd(f);
  d2sT2d_numerical = (Medium.EoS.dsTd(f_Tplus_dconst)-Medium.EoS.dsTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2s/dT2)@d=const analytical1= " + String(d2sT2d_analytical1));
  Modelica.Utilities.Streams.print("  (d2s/dT2)@d=const analytical2= " + String(d2sT2d_analytical2));
  Modelica.Utilities.Streams.print("  (d2s/dT2)@d=const   numerical= " + String(d2sT2d_numerical));
  // check (d2s/dT dd)
  d2sTd_analytical1 = Medium.EoS.d2sTd(f);
  d2sTd_analytical2 = -1/state.d^2*Medium.EoS.d2pT2d(f);
  d2sTd_numerical = (Medium.EoS.dsTd(f_dplus_Tconst)-Medium.EoS.dsTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2s/dT dd) analytical1= " + String(d2sTd_analytical1));
  Modelica.Utilities.Streams.print("  (d2s/dT dd) analytical2= " + String(d2sTd_analytical2));
  Modelica.Utilities.Streams.print("  (d2s/dT dd)   numerical= " + String(d2sTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("internal Energy wrt. dT");
  // check (du/dd)@T=const
  dudT_analytical1 = Medium.EoS.dudT(f);
  dudT_analytical2 = -state.T/state.d^2*Medium.EoS.dpTd(f) + state.p/state.d^2;
  dudT_numerical = (dplus_Tconst.u-dminus_Tconst.u)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (du/dd)@T=const analytical1= " + String(dudT_analytical1));
  Modelica.Utilities.Streams.print("  (du/dd)@T=const analytical2= " + String(dudT_analytical2));
  Modelica.Utilities.Streams.print("  (du/dd)@T=const   numerical= " + String(dudT_numerical));
  // check (du/dT)@d=const
  duTd_analytical = Medium.EoS.duTd(f);
  duTd_numerical = (Tplus_dconst.u-Tminus_dconst.u)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (du/dT)@d=const analytical= " + String(duTd_analytical));
  Modelica.Utilities.Streams.print("  (du/dT)@d=const  numerical= " + String(duTd_numerical));
  // check (d2u/dd2)@T=const
  d2ud2T_analytical1 = Medium.EoS.d2ud2T(f);
  d2ud2T_analytical2 = -state.T/state.d^2*Medium.EoS.d2pTd(f) + 2*state.T/state.d^3*Medium.EoS.dpTd(f) + 1/state.d^2*Medium.EoS.dpdT(f) -2*state.p/state.d^3;
  d2ud2T_numerical = (Medium.EoS.dudT(f_dplus_Tconst)-Medium.EoS.dudT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2u/dd2)@T=const analytical1= " + String(d2ud2T_analytical1));
  Modelica.Utilities.Streams.print("  (d2u/dd2)@T=const analytical2= " + String(d2ud2T_analytical2));
  Modelica.Utilities.Streams.print("  (d2u/dd2)@T=const   numerical= " + String(d2ud2T_numerical));
  // check (d2u/dT2)@d=const
  d2uT2d_analytical1 = Medium.EoS.d2uT2d(f);
  d2uT2d_analytical2 = state.T*Medium.EoS.d2sT2d(f) + Medium.EoS.dsTd(f);
  d2uT2d_numerical = (Medium.EoS.duTd(f_Tplus_dconst)-Medium.EoS.duTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2u/dT2)@d=const analytical1= " + String(d2uT2d_analytical1));
  Modelica.Utilities.Streams.print("  (d2u/dT2)@d=const analytical2= " + String(d2uT2d_analytical2));
  Modelica.Utilities.Streams.print("  (d2u/dT2)@d=const   numerical= " + String(d2uT2d_numerical));
  // check (d2u/dT dd)
  d2uTd_analytical1 = Medium.EoS.d2uTd(f);
  d2uTd_analytical2 = -state.T/state.d^2*Medium.EoS.d2pT2d(f);
  d2uTd_numerical = (Medium.EoS.duTd(f_dplus_Tconst)-Medium.EoS.duTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2u/dT dd) analytical1= " + String(d2uTd_analytical1));
  Modelica.Utilities.Streams.print("  (d2u/dT dd) analytical2= " + String(d2uTd_analytical2));
  Modelica.Utilities.Streams.print("  (d2u/dT dd)   numerical= " + String(d2uTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Enthalpy wrt. dT");
  // check (dh/dd)@T=const
  dhdT_analytical1 = Medium.EoS.dhdT(f);
  dhdT_analytical2 = 1/state.d*Medium.EoS.dpdT(f) - state.T/state.d^2*Medium.EoS.dpTd(f);
  dhdT_numerical = (dplus_Tconst.h-dminus_Tconst.h)/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const analytical1= " + String(dhdT_analytical1));
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const analytical2= " + String(dhdT_analytical2));
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const   numerical= " + String(dhdT_numerical));
  // check (dh/dT)@d=const
  dhTd_analytical1 = Medium.EoS.dhTd(f);
  dhTd_analytical2 = Medium.EoS.duTd(f) + 1/state.d*Medium.EoS.dpTd(f);
  dhTd_numerical = (Tplus_dconst.h-Tminus_dconst.h)/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const analytical1= " + String(dhTd_analytical1));
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const analytical2= " + String(dhTd_analytical2));
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const   numerical= " + String(dhTd_numerical));
  // check (d2h/dd2)@T=const
  d2hd2T_analytical1 = Medium.EoS.d2hd2T(f);
  d2hd2T_analytical2 = -state.T/state.d^2*Medium.EoS.d2pTd(f) + 1/state.d*Medium.EoS.d2pd2T(f) - 1/state.d^2*Medium.EoS.dpdT(f) + 2*state.T/state.d^3*Medium.EoS.dpTd(f);
  d2hd2T_numerical = (Medium.EoS.dhdT(f_dplus_Tconst)-Medium.EoS.dhdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2h/dd2)@T=const analytical1= " + String(d2hd2T_analytical1));
  Modelica.Utilities.Streams.print("  (d2h/dd2)@T=const analytical2= " + String(d2hd2T_analytical2));
  Modelica.Utilities.Streams.print("  (d2h/dd2)@T=const   numerical= " + String(d2hd2T_numerical));
  // check (d2h/dT2)@d=const
  d2hT2d_analytical1 = Medium.EoS.d2hT2d(f);
  d2hT2d_analytical2 = Medium.EoS.d2uT2d(f) + 1/state.d*Medium.EoS.d2pT2d(f);
  d2hT2d_numerical = (Medium.EoS.dhTd(f_Tplus_dconst)-Medium.EoS.dhTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2h/dT2)@d=const analytical1= " + String(d2hT2d_analytical1));
  Modelica.Utilities.Streams.print("  (d2h/dT2)@d=const analytical2= " + String(d2hT2d_analytical2));
  Modelica.Utilities.Streams.print("  (d2h/dT2)@d=const   numerical= " + String(d2hT2d_numerical));
  // check (d2h/dT dd)
  d2hTd_analytical1 = Medium.EoS.d2hTd(f);
  d2hTd_analytical2 = 1/state.d*Medium.EoS.d2pTd(f) - 1/state.d^2*Medium.EoS.dpTd(f) - T/state.d^2*Medium.EoS.d2pT2d(f);
  d2hTd_numerical = (Medium.EoS.dhTd(f_dplus_Tconst)-Medium.EoS.dhTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2h/dT dd) analytical1= " + String(d2hTd_analytical1));
  Modelica.Utilities.Streams.print("  (d2h/dT dd) analytical2= " + String(d2hTd_analytical2));
  Modelica.Utilities.Streams.print("  (d2h/dT dd)   numerical= " + String(d2hTd_numerical));
  Modelica.Utilities.Streams.print("Enthalpy wrt. pT");
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
  Modelica.Utilities.Streams.print("Gibbs energy wrt. dT");
  // check (dg/dd)@T=const
  dgdT_analytical1 = Medium.EoS.dgdT(f);
  dgdT_analytical2 = 1/state.d*Medium.EoS.dpdT(f);
  dgdT_numerical = ((dplus_Tconst.h-dplus_Tconst.T*dplus_Tconst.s)-(dminus_Tconst.h-dminus_Tconst.T*dminus_Tconst.s))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (dg/dd)@T=const analytical1= " + String(dgdT_analytical1));
  Modelica.Utilities.Streams.print("  (dg/dd)@T=const analytical2= " + String(dgdT_analytical2));
  Modelica.Utilities.Streams.print("  (dg/dd)@T=const   numerical= " + String(dgdT_numerical));
  // check (dg/dT)@d=const
  dgTd_analytical1 = Medium.EoS.dgTd(f);
  dgTd_analytical2 = -state.s + 1/state.d*Medium.EoS.dpTd(f);
  dgTd_numerical = ((Tplus_dconst.h-Tplus_dconst.T*Tplus_dconst.s)-(Tminus_dconst.h-Tminus_dconst.T*Tminus_dconst.s))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (dg/dT)@d=const analytical1= " + String(dgTd_analytical1));
  Modelica.Utilities.Streams.print("  (dg/dT)@d=const analytical2= " + String(dgTd_analytical2));
  Modelica.Utilities.Streams.print("  (dg/dT)@d=const   numerical= " + String(dgTd_numerical));
  // check (d2g/dd2)@T=const
  d2gd2T_analytical1 = Medium.EoS.d2gd2T(f);
  d2gd2T_analytical2 = 1/state.d*Medium.EoS.d2pd2T(f) - 1/state.d^2*Medium.EoS.dpdT(f);
  d2gd2T_numerical = (Medium.EoS.dgdT(f_dplus_Tconst)-Medium.EoS.dgdT(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2g/dd2)@T=const analytical1= " + String(d2gd2T_analytical1));
  Modelica.Utilities.Streams.print("  (d2g/dd2)@T=const analytical2= " + String(d2gd2T_analytical2));
  Modelica.Utilities.Streams.print("  (d2g/dd2)@T=const   numerical= " + String(d2gd2T_numerical));
  // check (d2g/dT2)@d=const
  d2gT2d_analytical1 = Medium.EoS.d2gT2d(f);
  d2gT2d_analytical2 = -1/state.T*Medium.EoS.duTd(f) +1/state.d*Medium.EoS.d2pT2d(f);
  d2gT2d_numerical = (Medium.EoS.dgTd(f_Tplus_dconst)-Medium.EoS.dgTd(f_Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d2g/dT2)@d=const analytical1= " + String(d2gT2d_analytical1));
  Modelica.Utilities.Streams.print("  (d2g/dT2)@d=const analytical2= " + String(d2gT2d_analytical2));
  Modelica.Utilities.Streams.print("  (d2g/dT2)@d=const   numerical= " + String(d2gT2d_numerical));
  // check (d2g/dT dd)
  d2gTd_analytical1 = Medium.EoS.d2gTd(f);
  d2gTd_analytical2 = 1/state.d*Medium.EoS.d2pTd(f);
  d2gTd_numerical = (Medium.EoS.dgTd(f_dplus_Tconst)-Medium.EoS.dgTd(f_dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d2g/dT dd) analytical1= " + String(d2gTd_analytical1));
  Modelica.Utilities.Streams.print("  (d2g/dT dd) analytical2= " + String(d2gTd_analytical2));
  Modelica.Utilities.Streams.print("  (d2g/dT dd)   numerical= " + String(d2gTd_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Density wrt. pT");
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
                      +(Medium.EoS.d2pTd(f)*Medium.EoS.dpdT(f) - Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f))/(Medium.EoS.dpdT(f)^3) * Medium.EoS.dpTd(f);
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
  Modelica.Utilities.Streams.print("Density wrt. ph");
  // check (dd/dp)@h=const
  ddph_analytical = Medium.density_derp_h(state=state);
  ddph_numerical = (pplus_hconst.d-pminus_hconst.d)/(pplus_hconst.p-pminus_hconst.p);
  Modelica.Utilities.Streams.print("  (dd/dp)@h=const analytical= " + String(ddph_analytical));
  Modelica.Utilities.Streams.print("  (dd/dp)@h=const  numerical= " + String(ddph_numerical));
  // check (dd/dh)@p=const
  ddhp_analytical = Medium.density_derh_p(state=state);
  ddhp_numerical = (hplus_pconst.d-hminus_pconst.d)/(hplus_pconst.h-hminus_pconst.h);
  Modelica.Utilities.Streams.print("  (dd/dh)@p=const analytical= " + String(ddhp_analytical));
  Modelica.Utilities.Streams.print("  (dd/dh)@p=const  numerical= " + String(ddhp_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Temperature wrt. pd");
  // check (dT/dd)@p=const
  dTdp_analytical = -Medium.EoS.dpdT(f)/Medium.EoS.dpTd(f);
  dTdp_numerical = (dplus_pconst.T-dminus_pconst.T)/(dplus_pconst.d-dminus_pconst.d);
  Modelica.Utilities.Streams.print("  (dT/dd)@p=const analytical= " + String(dTdp_analytical));
  Modelica.Utilities.Streams.print("  (dT/dd)@p=const  numerical= " + String(dTdp_numerical));
  // check (dT/dp)@d=const
  dTpd_analytical = 1/Medium.EoS.dpTd(f);
  dTpd_numerical = (pplus_dconst.T-pminus_dconst.T)/(pplus_dconst.p-pminus_dconst.p);
  Modelica.Utilities.Streams.print("  (dT/dp)@d=const analytical= " + String(dTpd_analytical));
  Modelica.Utilities.Streams.print("  (dT/dp)@d=const  numerical= " + String(dTpd_numerical));
  // check (d2T/dd2)@p=const
  d2Td2p_analytical = -(Medium.EoS.d2pd2T(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pTd(f))/(Medium.EoS.dpTd(f)^2)
                      +(Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pT2d(f))/(Medium.EoS.dpTd(f)^3) * Medium.EoS.dpdT(f);
  d2Td2p_numerical = (-Medium.EoS.dpdT(f_dplus_pconst)/Medium.EoS.dpTd(f_dplus_pconst)+Medium.EoS.dpdT(f_dminus_pconst)/Medium.EoS.dpTd(f_dminus_pconst))/(dplus_pconst.d-dminus_pconst.d);
  Modelica.Utilities.Streams.print("  (d2T/dd2)@p=const analytical= " + String(d2Td2p_analytical));
  Modelica.Utilities.Streams.print("  (d2T/dd2)@p=const  numerical= " + String(d2Td2p_numerical));
  // check (d2T/dp2)@d=const
  d2Tp2d_analytical = -Medium.EoS.d2pT2d(f)/Medium.EoS.dpTd(f)^3;
  d2Tp2d_numerical = (1/Medium.EoS.dpTd(f_pplus_dconst)-1/Medium.EoS.dpTd(f_pminus_dconst))/(pplus_dconst.p-pminus_dconst.p);
  Modelica.Utilities.Streams.print("  (d2T/dp2)@d=const analytical= " + String(d2Tp2d_analytical));
  Modelica.Utilities.Streams.print("  (d2T/dp2)@d=const  numerical= " + String(d2Tp2d_numerical));
  // check (d2T/dp dd)
  d2Tpd_analytical = -(Medium.EoS.d2pTd(f)*Medium.EoS.dpTd(f) - Medium.EoS.dpdT(f)*Medium.EoS.d2pT2d(f))/(Medium.EoS.dpTd(f)^3);
  d2Tpd_numerical1 = (1/Medium.EoS.dpTd(f_dplus_pconst)-1/Medium.EoS.dpTd(f_dminus_pconst))/(dplus_pconst.d-dminus_pconst.d);
  d2Tpd_numerical2 = (-Medium.EoS.dpdT(f_pplus_dconst)/Medium.EoS.dpTd(f_pplus_dconst)+Medium.EoS.dpdT(f_pminus_dconst)/Medium.EoS.dpTd(f_pminus_dconst))/(pplus_dconst.p-pminus_dconst.p);
  Modelica.Utilities.Streams.print("  (d2T/dp dd) analytical= " + String(d2Tpd_analytical));
  Modelica.Utilities.Streams.print("  (d2T/dp dd) numerical1= " + String(d2Tpd_numerical1));
  Modelica.Utilities.Streams.print("  (d2T/dd dp) numerical2= " + String(d2Tpd_numerical2));
  Modelica.Utilities.Streams.print("Temperature wrt. ph");
  // check (dT/dp)@h=const
  dTph_analytical = Medium.jouleThomsonCoefficient(state=state);
  dTph_numerical = (pplus_hconst.T-pminus_hconst.T)/(pplus_hconst.p-pminus_hconst.p);
  Modelica.Utilities.Streams.print("  (dT/dp)@h=const analytical= " + String(dTph_analytical));
  Modelica.Utilities.Streams.print("  (dT/dp)@h=const  numerical= " + String(dTph_numerical));

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print("Second order derivatives wrt. mixed properties");
  // check (d cp/dT)@d=const
  dcpTd_analytical =  Medium.EoS.d2hT2d(f) + Medium.EoS.dhdT(f)*Medium.EoS.dpTd(f)*Medium.EoS.d2pTd(f)/Medium.EoS.dpdT(f)^2
                   - (Medium.EoS.d2hTd(f)*Medium.EoS.dpTd(f) + Medium.EoS.dhdT(f)*Medium.EoS.d2pT2d(f))/Medium.EoS.dpdT(f);
  //dcpTd_analytical2 =  Medium.EoS.d2hT2d(f) - Medium.EoS.d2hTd(f)*Medium.EoS.dpTd(f)/Medium.EoS.dpdT(f);
  dcpTd_numerical = (Medium.specificHeatCapacityCp(Tplus_dconst)-Medium.specificHeatCapacityCp(Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d cp/dT)@d=const analytical= " + String(dcpTd_analytical));
  //Modelica.Utilities.Streams.print("  (d cp/dT)@d=const analytical2= " + String(dcpTd_analytical2));
  Modelica.Utilities.Streams.print("  (d cp/dT)@d=const  numerical= " + String(dcpTd_numerical));
  // check (d cp/dd)@T=const
  dcpdT_analytical =  Medium.EoS.d2hTd(f) + Medium.EoS.dhdT(f)*Medium.EoS.dpTd(f)*Medium.EoS.d2pd2T(f)/Medium.EoS.dpdT(f)^2
                   - (Medium.EoS.d2hd2T(f)*Medium.EoS.dpTd(f) + Medium.EoS.dhdT(f)*Medium.EoS.d2pTd(f))/Medium.EoS.dpdT(f);
  dcpdT_numerical = (Medium.specificHeatCapacityCp(dplus_Tconst)-Medium.specificHeatCapacityCp(dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d cp/dd)@T=const analytical= " + String(dcpdT_analytical));
  Modelica.Utilities.Streams.print("  (d cp/dd)@T=const  numerical= " + String(dcpdT_numerical));
  // check (d cp/dp)@h=const
  dcpph_analytical =  (dcpTd_analytical*Medium.EoS.dhdT(f) - dcpdT_analytical*Medium.EoS.dhTd(f))/(Medium.EoS.dpTd(f)*Medium.EoS.dhdT(f) - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f));
  dcpph_numerical = (Medium.specificHeatCapacityCp(pplus_hconst)-Medium.specificHeatCapacityCp(pminus_hconst))/(pplus_hconst.p-pminus_hconst.p);
  Modelica.Utilities.Streams.print("  (d cp/dp)@h=const analytical= " + String(dcpph_analytical));
  Modelica.Utilities.Streams.print("  (d cp/dp)@h=const  numerical= " + String(dcpph_numerical));

  // check (d mu/dT)@d=const
  dmuTd_analytical1 = (- Medium.EoS.d2pT2d(f) + (Medium.EoS.d2pTd(f)*Medium.EoS.dhTd(f) + Medium.EoS.dpdT(f)*Medium.EoS.d2hT2d(f))/Medium.EoS.dhdT(f)
                       - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f)*Medium.EoS.d2hTd(f)/Medium.EoS.dhdT(f)^2)
                       /(Medium.EoS.dpTd(f)-Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f)/Medium.EoS.dhdT(f))^2;
  dmuTd_analytical2 = (- Medium.EoS.d2pT2d(f) + (Medium.EoS.d2pTd(f)*Medium.EoS.dhTd(f) + Medium.EoS.dpdT(f)*Medium.EoS.d2hT2d(f))/Medium.EoS.dhdT(f)
                       - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f)*Medium.EoS.d2hTd(f)/Medium.EoS.dhdT(f)^2)
                       * Medium.jouleThomsonCoefficient(state)^2;
  dmuTd_numerical = (Medium.jouleThomsonCoefficient(Tplus_dconst)-Medium.jouleThomsonCoefficient(Tminus_dconst))/(Tplus_dconst.T-Tminus_dconst.T);
  Modelica.Utilities.Streams.print("  (d mu/dT)@d=const analytical1= " + String(dmuTd_analytical1));
  Modelica.Utilities.Streams.print("  (d mu/dT)@d=const analytical2= " + String(dmuTd_analytical2));
  Modelica.Utilities.Streams.print("  (d mu/dT)@d=const   numerical= " + String(dmuTd_numerical));
  // check (d mu/dd)@T=const
  dmudT_analytical1 = (- Medium.EoS.d2pTd(f) + (Medium.EoS.d2pd2T(f)*Medium.EoS.dhTd(f) + Medium.EoS.dpdT(f)*Medium.EoS.d2hTd(f))/Medium.EoS.dhdT(f)
                       - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f)*Medium.EoS.d2hd2T(f)/Medium.EoS.dhdT(f)^2)
                       /(Medium.EoS.dpTd(f)-Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f)/Medium.EoS.dhdT(f))^2;
  dmudT_analytical2 = (- Medium.EoS.d2pTd(f) + (Medium.EoS.d2pd2T(f)*Medium.EoS.dhTd(f) + Medium.EoS.dpdT(f)*Medium.EoS.d2hTd(f))/Medium.EoS.dhdT(f)
                       - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f)*Medium.EoS.d2hd2T(f)/Medium.EoS.dhdT(f)^2)
                       * Medium.jouleThomsonCoefficient(state)^2;
  dmudT_numerical = (Medium.jouleThomsonCoefficient(dplus_Tconst)-Medium.jouleThomsonCoefficient(dminus_Tconst))/(dplus_Tconst.d-dminus_Tconst.d);
  Modelica.Utilities.Streams.print("  (d mu/dd)@T=const analytical1= " + String(dmudT_analytical1));
  Modelica.Utilities.Streams.print("  (d mu/dd)@T=const analytical2= " + String(dmudT_analytical2));
  Modelica.Utilities.Streams.print("  (d mu/dd)@T=const   numerical= " + String(dmudT_numerical));
  // check (d mu/dp)@h=const
  dmuph_analytical1 = (dmuTd_analytical1*Medium.EoS.dhdT(f) - dmudT_analytical1*Medium.EoS.dhTd(f))/(Medium.EoS.dpTd(f)*Medium.EoS.dhdT(f) - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f));
  dmuph_analytical2 = (dmuTd_analytical2*Medium.EoS.dhdT(f) - dmudT_analytical2*Medium.EoS.dhTd(f))/(Medium.EoS.dpTd(f)*Medium.EoS.dhdT(f) - Medium.EoS.dpdT(f)*Medium.EoS.dhTd(f));
  dmuph_numerical = (Medium.jouleThomsonCoefficient(pplus_hconst)-Medium.jouleThomsonCoefficient(pminus_hconst))/(pplus_hconst.p-pminus_hconst.p);
  Modelica.Utilities.Streams.print("  (d mu/dp)@h=const analytical1= " + String(dmuph_analytical1));
  Modelica.Utilities.Streams.print("  (d mu/dp)@h=const analytical2= " + String(dmuph_analytical2));
  Modelica.Utilities.Streams.print("  (d mu/dp)@h=const   numerical= " + String(dmuph_numerical));

annotation (experiment(NumberOfIntervals=1));
end Derivatives_SinglePhase;
