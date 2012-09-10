within HelmholtzMedia.Examples.Validation;
model Derivatives_TwoPhase
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Butane;
  // choose d and T which will result in two-phase
  parameter Medium.Density d=228;
  parameter Medium.Temperature T=220;
  Medium.ThermodynamicState state=Medium.setState_dTX(d=d, T=T);
  Medium.SaturationProperties sat=Medium.setSat_T(T=T);

// Vapour mass fraction derivatives
  Real dxTv_numerical;
  Real dxTv_analytical1;
  Real dxTv_analytical2;
  Real dxTv_analytical3;
  Real dxTv_analytical4;
  Real dxph_numerical;
  Real dxph_analytical1;
  Real dxph_analytical2;
  /*Real dxTs_numerical;
  Real dxTs_analytical;
  Real dxTs_analytical2;
  Real dxTh_numerical;
  Real dxTh_analytical;
  Real dxTh_analytical2; */
// Entropy derivatives
  Medium.Types.DerEntropyByTemperature dsTd_numerical;
  Medium.Types.DerEntropyByTemperature dsTd_analytical;
// Energy derivatives
  Medium.Types.DerEnergyByTemperature duTd_numerical;
  Medium.Types.DerEnergyByTemperature duTd_analytical1;
  Medium.Types.DerEnergyByTemperature duTd_analytical2;
  Medium.Types.DerEnergyByTemperature duTd_analytical3;
// Density derivatives
  Medium.DerDensityByEnthalpy ddhp_numerical;
  Medium.DerDensityByEnthalpy ddhp_analytical1;
  Medium.DerDensityByEnthalpy ddhp_analytical2;
  Medium.DerDensityByEnthalpy ddhp_analytical3;
  Medium.DerDensityByPressure ddph_numerical;
  Medium.DerDensityByPressure ddph_analytical1;
  Medium.DerDensityByPressure ddph_analytical2;
// Enthalpy derivatives
  Medium.Types.DerEnthalpyByDensity dhdT_numerical;
  Medium.Types.DerEnthalpyByDensity dhdT_analytical;
  Medium.Types.DerEnthalpyByTemperature dhTd_numerical;
  Medium.Types.DerEnthalpyByTemperature dhTd_analytical;

protected
  Medium.Types.DerPressureByTemperature dpT = (sat.vap.s-sat.liq.s)/(1.0/sat.vap.d-1.0/sat.liq.d);
  Medium.Types.DerTemperatureByPressure dTp = (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);
  Medium.EoS.HelmholtzDerivs fl=Medium.EoS.setHelmholtzDerivsSecond(T=T, d=sat.liq.d, phase=1);
  Medium.EoS.HelmholtzDerivs fv=Medium.EoS.setHelmholtzDerivsSecond(T=T, d=sat.vap.d, phase=1);
  Medium.MassFraction x=Medium.vapourQuality(state=state);
  Medium.ThermodynamicState d_plus=Medium.setState_dTX(d=d*1.0001, T=T);
  Medium.ThermodynamicState d_minus=Medium.setState_dTX(d=d*0.9999, T=T);
  Medium.ThermodynamicState T_plus=Medium.setState_dTX(d=d, T=T*1.0001);
  Medium.ThermodynamicState T_minus=Medium.setState_dTX(d=d, T=T*0.9999);
  Medium.ThermodynamicState h_plus=Medium.setState_phX(p=state.p, h=state.h+abs(0.0001*state.h));
  Medium.ThermodynamicState h_minus=Medium.setState_phX(p=state.p, h=state.h-abs(0.0001*state.h));
  Medium.ThermodynamicState p_plus=Medium.setState_phX(p=state.p*1.0001, h=state.h);
  Medium.ThermodynamicState p_minus=Medium.setState_phX(p=state.p*0.9999, h=state.h);
// Entropy derivatives along saturation
  Medium.Types.DerEntropyByTemperature dsTd_liq = fl.R/T*(-fl.tau^2*(fl.itt+fl.rtt));
  Medium.Types.DerEntropyByDensity dsdT_liq = fl.R/sat.liq.d*(-(1+fl.delta*fl.rd)+(0+fl.tau*fl.delta*fl.rtd));
  Medium.Types.DerEntropyByTemperature dsTp_liq = dsTd_liq-dsdT_liq*Medium.pressure_derT_d(state=sat.liq)/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEntropyByPressure dspT_liq = dsdT_liq/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEntropyByTemperature dsT_liq = dsTp_liq+dspT_liq*dpT;
  Medium.Types.DerEntropyByPressure dsp_liq = dspT_liq+dsTp_liq*Medium.saturationTemperature_derp(p=state.p);
  Medium.Types.DerEntropyByTemperature dsTd_vap = fv.R/T*(-fv.tau^2*(fv.itt+fv.rtt));
  Medium.Types.DerEntropyByDensity dsdT_vap = fv.R/sat.vap.d*(-(1+fv.delta*fv.rd)+(0+fv.tau*fv.delta*fv.rtd));
  Medium.Types.DerEntropyByTemperature dsTp_vap = dsTd_vap-dsdT_vap*Medium.pressure_derT_d(state=sat.vap)/Medium.pressure_derd_T(state=sat.vap);
  Medium.Types.DerEntropyByPressure dspT_vap = dsdT_vap/Medium.pressure_derd_T(state=sat.vap);
  Medium.Types.DerEntropyByTemperature dsT_vap = dsTp_vap+dspT_vap*dpT;
  Medium.Types.DerEntropyByPressure dsp_vap = dspT_vap+dsTp_vap*dTp;
// Internal energy derivatives along saturation line
  Medium.Types.DerEnergyByDensity dudT_liq = fl.R*T/sat.liq.d*fl.tau*fl.delta*fl.rtd;
  Medium.Types.DerEnergyByTemperature duTd_liq = Medium.specificHeatCapacityCv(state=sat.liq);
  Medium.Types.DerEnergyByTemperature duTp_liq = duTd_liq-dudT_liq*Medium.pressure_derT_d(state=sat.liq)/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEnergyByPressure dupT_liq = dudT_liq/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEnergyByTemperature duT_liq = duTp_liq+dupT_liq*dpT;
  Medium.Types.DerEnergyByDensity dudT_vap = fv.R*T/sat.vap.d*fv.tau*fv.delta*fv.rtd;
  Medium.Types.DerEnergyByTemperature duTd_vap = Medium.specificHeatCapacityCv(state=sat.vap);
  Medium.Types.DerEnergyByTemperature duTp_vap = duTd_vap-dudT_vap*Medium.pressure_derT_d(state=sat.vap)/Medium.pressure_derd_T(state=sat.vap);
  Medium.Types.DerEnergyByPressure dupT_vap = dudT_vap/Medium.pressure_derd_T(state=sat.vap);
  Medium.Types.DerEnergyByTemperature duT_vap = duTp_vap+dupT_vap*dpT;
// Density derivatives along saturation line
  Medium.DerDensityByTemperature ddT_liq = Medium.density_derT_p(state=sat.liq) +Medium.density_derp_T(state=sat.liq)*dpT;
  Medium.DerDensityByTemperature ddT_vap = Medium.density_derT_p(state=sat.vap) +Medium.density_derp_T(state=sat.vap)*dpT;
  Medium.DerDensityByPressure ddp_liq = Medium.density_derp_T(state=sat.liq) +Medium.density_derT_p(state=sat.liq)*dTp;
  Medium.DerDensityByPressure ddp_vap = Medium.density_derp_T(state=sat.vap) +Medium.density_derT_p(state=sat.vap)*dTp;
// Specific Volume derivatives along saturatin line
  Real dvT_liq = -1/sat.liq.d^2 * ddT_liq;
  Real dvT_vap = -1/sat.vap.d^2 * ddT_vap;
  Real dvp_liq = -1/sat.liq.d^2 * ddp_liq;
  Real dvp_vap = -1/sat.vap.d^2 * ddp_vap;
// Enthalpy derivatives along saturation line
  Medium.Types.DerEnthalpyByTemperature dhT_liq = Medium.specificHeatCapacityCp(state=sat.liq) + Medium.isothermalThrottlingCoefficient(state=sat.liq)*dpT;
  Medium.Types.DerEnthalpyByTemperature dhT_vap = Medium.specificHeatCapacityCp(state=sat.vap) + Medium.isothermalThrottlingCoefficient(state=sat.vap)*dpT;
  Medium.Types.DerEnthalpyByPressure dhp_liq = Medium.dBubbleEnthalpy_dPressure(sat=sat);
  Medium.Types.DerEnthalpyByPressure dhp_vap = Medium.dDewEnthalpy_dPressure(sat=sat);

equation
  assert(state.phase == 2, "state not in two-phase region");

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Vapour mass fraction");
  dxTv_numerical = (Medium.vapourQuality(T_plus) - Medium.vapourQuality(T_minus))/(T_plus.T - T_minus.T);
  dxTv_analytical1 = (-dvT_liq*(1/sat.vap.d - 1/sat.liq.d) - (1/state.d - 1/sat.liq.d)*(dvT_vap - dvT_liq))/(1/sat.vap.d - 1/sat.liq.d)^2;
  dxTv_analytical2 = (+1/sat.liq.d^2*ddT_liq*(1/sat.vap.d - 1/sat.liq.d) - (1/state.d - 1/sat.liq.d)*(-1/sat.vap.d^2*ddT_vap + 1/sat.liq.d^2*ddT_liq))/(1/sat.vap.d - 1/sat.liq.d)^2;
  dxTv_analytical3 = (x*dvT_vap + (1 - x)*dvT_liq)/(1/sat.liq.d - 1/sat.vap.d);
  dxTv_analytical4 = (x*(-1)/sat.vap.d^2*ddT_vap + (1 - x)*(-1)/sat.liq.d^2*ddT_liq)/(1/sat.liq.d - 1/sat.vap.d);
  Modelica.Utilities.Streams.print("  (dx/dT)@v=const   numerical= " + String(dxTv_numerical));
  Modelica.Utilities.Streams.print("  (dx/dT)@v=const analytical1= " + String(dxTv_analytical1));
  Modelica.Utilities.Streams.print("  (dx/dT)@v=const analytical2= " + String(dxTv_analytical2));
  Modelica.Utilities.Streams.print("  (dx/dT)@v=const analytical3= " + String(dxTv_analytical3));
  Modelica.Utilities.Streams.print("  (dx/dT)@v=const analytical4= " + String(dxTv_analytical4));

  Modelica.Utilities.Streams.print(" ");
  dxph_numerical = (Medium.vapourQuality(p_plus) - Medium.vapourQuality(p_minus))/(p_plus.p - p_minus.p);
  dxph_analytical1 = (-dhp_liq*(sat.vap.h - sat.liq.h) - (state.h - sat.liq.h)*(dhp_vap - dhp_liq))/(sat.vap.h - sat.liq.h)^2;
  dxph_analytical2 = (x*dhp_vap + (1 - x)*dhp_liq)/(sat.liq.h - sat.vap.h);
  Modelica.Utilities.Streams.print("  (dx/dp)@h=const   numerical= " + String(dxph_numerical));
  Modelica.Utilities.Streams.print("  (dx/dp)@h=const analytical1= " + String(dxph_analytical1));
  Modelica.Utilities.Streams.print("  (dx/dp)@h=const analytical2= " + String(dxph_analytical2));

  // dxTh_analytical1 = (-dhT_liq *(sat.vap.h-sat.liq.h) -(state.h-sat.liq.h)*(dhT_vap -dhT_liq) ) / (sat.vap.h-sat.liq.h)^2;
  // dxTs_analytical1 = (-dsT_liq *(sat.vap.s-sat.liq.s) -(state.s-sat.liq.s)*(dsT_vap -dsT_liq) ) / (sat.vap.s-sat.liq.s)^2;
  // dxTh_analytical2 = (x*dhT_vap+(1-x)*dhT_liq) / (sat.liq.h-sat.vap.h);
  // dxTs_analytical2 = (x*dsT_vap+(1-x)*dsT_liq) / (sat.liq.s-sat.vap.s);
  // dxps_analytical1 = (-dsp_liq*(sat.vap.s-sat.liq.s)-(state.s-sat.liq.s)*(dsp_vap-dsp_liq)) / (sat.vap.s-sat.liq.s)^2;

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Entropy");
  // check (ds/dT)@d=const
  dsTd_numerical = (T_plus.s - T_minus.s)/(T_plus.T - T_minus.T);
  dsTd_analytical = (dsT_liq + x*(dsT_vap - dsT_liq) + dxTv_analytical1*(sat.vap.s-sat.liq.s));
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const  numerical= " + String(dsTd_numerical));
  Modelica.Utilities.Streams.print("  (ds/dT)@d=const analytical= " + String(dsTd_analytical));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Internal energy");
  // check (du/dd)@T=const
  // dudT_analytical = f.R*T/d*f.tau*f.delta*f.rtd;
  // dudT_numerical = (d_plus.u-d_minus.u)/(d_plus.d-d_minus.d);
  // Modelica.Utilities.Streams.print("(du/dd)@T=const analytical= " + String(dudT_analytical));
  // Modelica.Utilities.Streams.print("(du/dd)@T=const  numerical= " + String(dudT_numerical));
  // check (du/dT)@d=const
  duTd_numerical = (T_plus.u - T_minus.u)/(T_plus.T - T_minus.T);
  duTd_analytical1 = Medium.specificHeatCapacityCv(state=state);
  duTd_analytical2 = (duT_liq + x*(duT_vap - duT_liq) + dxTv_analytical1*(sat.vap.u-sat.liq.u));
  duTd_analytical3 = T*dsTd_analytical;
  Modelica.Utilities.Streams.print("  (du/dT)@d=const   numerical= " + String(duTd_numerical));
  Modelica.Utilities.Streams.print("  (du/dT)@d=const analytical1= " + String(duTd_analytical1));
  Modelica.Utilities.Streams.print("  (du/dT)@d=const analytical2= " + String(duTd_analytical2));
  Modelica.Utilities.Streams.print("  (du/dT)@d=const analytical3= " + String(duTd_analytical3));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Density");
  // check (dd/dh)@p=const
  ddhp_numerical = (h_plus.d - h_minus.d)/(h_plus.h - h_minus.h);
  ddhp_analytical1 = Medium.density_derh_p(state=state);
  ddhp_analytical2 = -state.d^2*(1/sat.liq.d - 1/sat.vap.d)/(sat.liq.h - sat.vap.h);
  ddhp_analytical3 = -state.d^2/T*dTp;
  Modelica.Utilities.Streams.print("  (dd/dh)@p=const   numerical= " + String(ddhp_numerical));
  Modelica.Utilities.Streams.print("  (dd/dh)@p=const analytical1= " + String(ddhp_analytical1));
  Modelica.Utilities.Streams.print("  (dd/dh)@p=const analytical2= " + String(ddhp_analytical2));
  Modelica.Utilities.Streams.print("  (dd/dh)@p=const analytical3= " + String(ddhp_analytical3));
  // check (dd/dp)@h=const
  ddph_numerical = (p_plus.d - p_minus.d)/(p_plus.p - p_minus.p);
  ddph_analytical1 = Medium.density_derp_h(state=state);
  ddph_analytical2 = -state.d^2*(dvp_liq + x*(dvp_vap - dvp_liq) + dxph_analytical1*(1/sat.vap.d - 1/sat.liq.d));
  Modelica.Utilities.Streams.print("  (dd/dp)@h=const   numerical= " + String(ddph_numerical));
  Modelica.Utilities.Streams.print("  (dd/dp)@h=const analytical1= " + String(ddph_analytical1));
  Modelica.Utilities.Streams.print("  (dd/dp)@h=const analytical2= " + String(ddph_analytical2));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Enthalpy");
  // check (dh/dd)@T=const
  dhdT_numerical = (d_plus.h - d_minus.h)/(d_plus.d - d_minus.d);
  dhdT_analytical = Medium.specificEnthalpy_derd_T(state=state);
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const  numerical= " + String(dhdT_numerical));
  Modelica.Utilities.Streams.print("  (dh/dd)@T=const analytical= " + String(dhdT_analytical));
  // check (dh/dT)@d=const
  dhTd_numerical = (T_plus.h - T_minus.h)/(T_plus.T - T_minus.T);
  dhTd_analytical = Medium.specificHeatCapacityCv(state=state) + 1/state.d*dpT;
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const  numerical= " + String(dhTd_numerical));
  Modelica.Utilities.Streams.print("  (dh/dT)@d=const analytical= " + String(dhTd_analytical));

annotation (experiment(NumberOfIntervals=1));
end Derivatives_TwoPhase;
