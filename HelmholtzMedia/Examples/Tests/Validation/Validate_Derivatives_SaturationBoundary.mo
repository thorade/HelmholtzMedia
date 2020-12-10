within HelmholtzMedia.Examples.Tests.Validation;
model Validate_Derivatives_SaturationBoundary
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Butane;
  parameter Medium.Temperature T=420;  // T has to be below T_crit

// Density derivatives
  Medium.DerDensityByTemperature ddT_liq_numerical;
  Medium.DerDensityByTemperature ddT_liq_analytical;
  Medium.DerDensityByTemperature ddT_vap_numerical;
  Medium.DerDensityByTemperature ddT_vap_analytical;
  Medium.DerDensityByPressure ddp_liq_numerical;
  Medium.DerDensityByPressure ddp_liq_analytical1;
  Medium.DerDensityByPressure ddp_liq_analytical2;
  Medium.DerDensityByPressure ddp_vap_numerical;
  Medium.DerDensityByPressure ddp_vap_analytical1;
  Medium.DerDensityByPressure ddp_vap_analytical2;
// Enthalpy derivatives
  Medium.Types.DerEnthalpyByTemperature dhT_liq_analytical;
  Medium.Types.DerEnthalpyByTemperature dhT_liq_numerical;
  Medium.Types.DerEnthalpyByTemperature dhT_vap_analytical;
  Medium.Types.DerEnthalpyByTemperature dhT_vap_numerical;
  Medium.Types.DerEnthalpyByPressure dhp_liq_analytical;
  Medium.Types.DerEnthalpyByPressure dhp_liq_numerical;
  Medium.Types.DerEnthalpyByPressure dhp_vap_analytical;
  Medium.Types.DerEnthalpyByPressure dhp_vap_numerical;
// Entropy derivatives
  Medium.Types.DerEntropyByTemperature dsT_liq_analytical;
  Medium.Types.DerEntropyByTemperature dsT_liq_numerical;
  Medium.Types.DerEntropyByTemperature dsT_vap_analytical;
  Medium.Types.DerEntropyByTemperature dsT_vap_numerical;
  Medium.Types.DerEntropyByPressure dsp_liq_analytical;
  Medium.Types.DerEntropyByPressure dsp_liq_numerical;
  Medium.Types.DerEntropyByPressure dsp_vap_analytical;
  Medium.Types.DerEntropyByPressure dsp_vap_numerical;
//saturated heat capacity
  Medium.Types.DerEnthalpyByTemperature c_sigma_liq_Span2000;
  Medium.Types.DerEnthalpyByTemperature c_sigma_liq_analytical;
  Medium.Types.DerEnthalpyByTemperature c_sigma_liq_numerical;
  Medium.Types.DerEnthalpyByTemperature c_sigma_vap_Span2000;
  Medium.Types.DerEnthalpyByTemperature c_sigma_vap_analytical;
  Medium.Types.DerEnthalpyByTemperature c_sigma_vap_numerical;
// Energy derivatives
  Medium.Types.DerEnergyByTemperature duT_liq_numerical;
  Medium.Types.DerEnergyByTemperature duT_liq_analytical;
  Medium.Types.DerEnergyByTemperature duT_vap_numerical;
  Medium.Types.DerEnergyByTemperature duT_vap_analytical;

protected
  Medium.SaturationProperties sat=Medium.setSat_T(T=T);
  Medium.HelmholtzDerivs fl=Medium.setHelmholtzDerivs(T=T, d=sat.liq.d, phase=1);
  Medium.HelmholtzDerivs fv=Medium.setHelmholtzDerivs(T=T, d=sat.vap.d, phase=1);
  Medium.SaturationProperties sat_Tplus = Medium.setSat_T(T=1.0001*T);
  Medium.SaturationProperties sat_Tminus= Medium.setSat_T(T=0.9999*T);
  Medium.SaturationProperties sat_pplus = Medium.setSat_p(p=1.0001*sat.psat);
  Medium.SaturationProperties sat_pminus= Medium.setSat_p(p=0.9999*sat.psat);
// Entropy single phase derivatives
  Medium.Types.DerEntropyByDensity dsdT_liq = fl.R_s/sat.liq.d*(-(1+fl.delta*fl.rd)+(0+fl.tau*fl.delta*fl.rtd));
  Medium.Types.DerEntropyByTemperature dsTd_liq = fl.R_s/T*(-fl.tau^2*(fl.itt+fl.rtt));
  Medium.Types.DerEntropyByTemperature dsTp_liq = dsTd_liq-dsdT_liq*Medium.pressure_derT_d(state=sat.liq)/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEntropyByPressure dspT_liq = dsdT_liq/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEntropyByDensity dsdT_vap = fv.R_s/sat.vap.d*(-(1+fv.delta*fv.rd)+(0+fv.tau*fv.delta*fv.rtd));
  Medium.Types.DerEntropyByTemperature dsTd_vap = fv.R_s/T*(-fv.tau^2*(fv.itt+fv.rtt));
  Medium.Types.DerEntropyByTemperature dsTp_vap = dsTd_vap-dsdT_vap*Medium.pressure_derT_d(state=sat.vap)/Medium.pressure_derd_T(state=sat.vap);
  Medium.Types.DerEntropyByPressure dspT_vap = dsdT_vap/Medium.pressure_derd_T(state=sat.vap);
// Internal energy single phase derivatives
  Medium.Types.DerEnergyByDensity dudT_liq = fl.R_s*T/sat.liq.d*fl.tau*fl.delta*fl.rtd;
  Medium.Types.DerEnergyByTemperature duTd_liq = Medium.specificHeatCapacityCv(state=sat.liq);
  Medium.Types.DerEnergyByTemperature duTp_liq = duTd_liq-dudT_liq*Medium.pressure_derT_d(state=sat.liq)/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEnergyByPressure dupT_liq = dudT_liq/Medium.pressure_derd_T(state=sat.liq);
  Medium.Types.DerEnergyByDensity dudT_vap = fv.R_s*T/sat.vap.d*fv.tau*fv.delta*fv.rtd;
  Medium.Types.DerEnergyByTemperature duTd_vap = Medium.specificHeatCapacityCv(state=sat.vap);
  Medium.Types.DerEnergyByTemperature duTp_vap = duTd_vap-dudT_vap*Medium.pressure_derT_d(state=sat.vap)/Medium.pressure_derd_T(state=sat.vap);
  Medium.Types.DerEnergyByPressure dupT_vap = dudT_vap/Medium.pressure_derd_T(state=sat.vap);

equation
  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  Modelica.Utilities.Streams.print(" ");

  Modelica.Utilities.Streams.print("Density");
  // check (dd/dT)@liq
  ddT_liq_numerical = (sat_Tplus.liq.d - sat_Tminus.liq.d)/(sat_Tplus.liq.T - sat_Tminus.liq.T);
  ddT_liq_analytical = Medium.density_derT_p(state=sat.liq) +Medium.density_derp_T(state=sat.liq)*Medium.saturationPressure_derT(T=sat.Tsat, sat=sat);
  Modelica.Utilities.Streams.print("  (dd/dT)@liq  numerical= " + String(ddT_liq_numerical));
  Modelica.Utilities.Streams.print("  (dd/dT)@liq analytical= " + String(ddT_liq_analytical));
  // check (dd/dT)@vap
  ddT_vap_numerical = (sat_Tplus.vap.d - sat_Tminus.vap.d)/(sat_Tplus.vap.T - sat_Tminus.vap.T);
  ddT_vap_analytical = Medium.density_derT_p(state=sat.vap) +Medium.density_derp_T(state=sat.vap)*Medium.saturationPressure_derT(T=sat.Tsat, sat=sat);
  Modelica.Utilities.Streams.print("  (dd/dT)@vap   numerical= " + String(ddT_vap_numerical));
  Modelica.Utilities.Streams.print("  (dd/dT)@vap analytical= " + String(ddT_vap_analytical));
  // check (dd/dp)@liq
  ddp_liq_numerical = (sat_pplus.liq.d - sat_pminus.liq.d)/(sat_pplus.liq.p - sat_pminus.liq.p);
  ddp_liq_analytical1 = Medium.density_derp_T(state=sat.liq) +Medium.density_derT_p(state=sat.liq)*Medium.saturationTemperature_derp(p=sat.psat, sat=sat);
  ddp_liq_analytical2 = Medium.dBubbleDensity_dPressure(sat=sat);
  Modelica.Utilities.Streams.print("  (dd/dp)@liq   numerical= " + String(ddp_liq_numerical));
  Modelica.Utilities.Streams.print("  (dd/dp)@liq analytical1= " + String(ddp_liq_analytical1));
  Modelica.Utilities.Streams.print("  (dd/dp)@liq analytical2= " + String(ddp_liq_analytical2));
  // check (dd/dp)@vap
  ddp_vap_numerical = (sat_pplus.vap.d - sat_pminus.vap.d)/(sat_pplus.vap.p - sat_pminus.vap.p);
  ddp_vap_analytical1 = Medium.density_derp_T(state=sat.vap) +Medium.density_derT_p(state=sat.vap)*Medium.saturationTemperature_derp(p=sat.psat, sat=sat);
  ddp_vap_analytical2 = Medium.dDewDensity_dPressure(sat=sat);
  Modelica.Utilities.Streams.print("  (dd/dp)@vap   numerical= " + String(ddp_vap_numerical));
  Modelica.Utilities.Streams.print("  (dd/dp)@vap analytical1= " + String(ddp_vap_analytical1));
  Modelica.Utilities.Streams.print("  (dd/dp)@vap analytical2= " + String(ddp_vap_analytical2));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Enthalpy");
  // check (dh/dT)@liq
  dhT_liq_numerical = (sat_Tplus.liq.h-sat_Tminus.liq.h)/(sat_Tplus.liq.T-sat_Tminus.liq.T);
  dhT_liq_analytical = Medium.specificHeatCapacityCp(state=sat.liq) +Medium.isothermalThrottlingCoefficient(state=sat.liq)*Medium.saturationPressure_derT(T=sat.Tsat, sat=sat);
  Modelica.Utilities.Streams.print("  (dh/dT)@liq  numerical= " + String(dhT_liq_numerical));
  Modelica.Utilities.Streams.print("  (dh/dT)@liq analytical= " + String(dhT_liq_analytical));
  // check (dh/dT)@vap
  dhT_vap_numerical = (sat_Tplus.vap.h-sat_Tminus.vap.h)/(sat_Tplus.vap.T-sat_Tminus.vap.T);
  dhT_vap_analytical = Medium.specificHeatCapacityCp(state=sat.vap) +Medium.isothermalThrottlingCoefficient(state=sat.vap)*Medium.saturationPressure_derT(T=sat.Tsat, sat=sat);
  Modelica.Utilities.Streams.print("  (dh/dT)@vap  numerical= " + String(dhT_vap_numerical));
  Modelica.Utilities.Streams.print("  (dh/dT)@vap analytical= " + String(dhT_vap_analytical));
  // check (dh/dp)@liq
  dhp_liq_numerical = (sat_pplus.liq.h-sat_pminus.liq.h)/(sat_pplus.liq.p-sat_pminus.liq.p);
  dhp_liq_analytical = Medium.dBubbleEnthalpy_dPressure(sat=sat);
  Modelica.Utilities.Streams.print("  (dh/dp)@liq  numerical= " + String(dhp_liq_numerical));
  Modelica.Utilities.Streams.print("  (dh/dp)@liq analytical= " + String(dhp_liq_analytical));
  // check (dh/dp)@vap
  dhp_vap_numerical = (sat_pplus.vap.h-sat_pminus.vap.h)/(sat_pplus.vap.p-sat_pminus.vap.p);
  dhp_vap_analytical = Medium.dDewEnthalpy_dPressure(sat=sat);
  Modelica.Utilities.Streams.print("  (dh/dp)@vap  numerical= " + String(dhp_vap_numerical));
  Modelica.Utilities.Streams.print("  (dh/dp)@vap analytical= " + String(dhp_vap_analytical));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Entropy");
  // check (ds/dT)@liq
  dsT_liq_numerical = (sat_Tplus.liq.s-sat_Tminus.liq.s)/(sat_Tplus.liq.T-sat_Tminus.liq.T);
  dsT_liq_analytical = dsTp_liq+dspT_liq*Medium.saturationPressure_derT(T=T);
  Modelica.Utilities.Streams.print("  (ds/dT)@liq  numerical= " + String(dsT_liq_numerical));
  Modelica.Utilities.Streams.print("  (ds/dT)@liq analytical= " + String(dsT_liq_analytical));
  // check (ds/dT)@vap
  dsT_vap_numerical = (sat_Tplus.vap.s-sat_Tminus.vap.s)/(sat_Tplus.vap.T-sat_Tminus.vap.T);
  dsT_vap_analytical = dsTp_vap+dspT_vap*Medium.saturationPressure_derT(T=T);
  Modelica.Utilities.Streams.print("  (ds/dT)@vap  numerical= " + String(dsT_vap_numerical));
  Modelica.Utilities.Streams.print("  (ds/dT)@vap analytical= " + String(dsT_vap_analytical));
  // check (ds/dp)@liq
  dsp_liq_numerical = (sat_pplus.liq.s-sat_pminus.liq.s)/(sat_pplus.liq.p-sat_pminus.liq.p);
  dsp_liq_analytical = dspT_liq+dsTp_liq*Medium.saturationTemperature_derp(p=sat.psat);
  Modelica.Utilities.Streams.print("  (ds/dp)@liq  numerical= " + String(dsp_liq_numerical));
  Modelica.Utilities.Streams.print("  (ds/dp)@liq analytical= " + String(dsp_liq_analytical));
  // check (ds/dp)@vap
  dsp_vap_numerical = (sat_pplus.vap.s-sat_pminus.vap.s)/(sat_pplus.vap.p-sat_pminus.vap.p);
  dsp_vap_analytical = dspT_vap+dsTp_vap*Medium.saturationTemperature_derp(p=sat.psat);
  Modelica.Utilities.Streams.print("  (ds/dp)@vap  numerical= " + String(dsp_vap_numerical));
  Modelica.Utilities.Streams.print("  (ds/dp)@vap analytical= " + String(dsp_vap_analytical));

  Modelica.Utilities.Streams.print(" ");
   Modelica.Utilities.Streams.print("saturated heat capacity");
  // check (c_sigma)@liq
  c_sigma_liq_numerical = T*dsT_liq_numerical;
  c_sigma_liq_analytical = T*dsT_liq_analytical;
  c_sigma_liq_Span2000 = Medium.specificHeatCapacityCv(state=sat.liq) - T*Medium.pressure_derT_d(state=sat.liq)/(sat.liq.d^2)*ddT_liq_analytical;
  Modelica.Utilities.Streams.print("  (c_sigma)@liq numerical= " + String(c_sigma_liq_numerical));
  Modelica.Utilities.Streams.print("  (c_sigma)@liq analytical= " + String(c_sigma_liq_analytical));
  Modelica.Utilities.Streams.print("  (c_sigma)@liq Span= " + String(c_sigma_liq_Span2000));
  // check (c_sigma)@vap
  c_sigma_vap_numerical = T*dsT_vap_numerical;
  c_sigma_vap_analytical = T*dsT_vap_analytical;
  c_sigma_vap_Span2000 = Medium.specificHeatCapacityCv(state=sat.vap) - T*Medium.pressure_derT_d(state=sat.vap)/(sat.vap.d^2)*ddT_vap_analytical;
  Modelica.Utilities.Streams.print("  (c_sigma)@vap numerical= " + String(c_sigma_vap_numerical));
  Modelica.Utilities.Streams.print("  (c_sigma)@vap analytical= " + String(c_sigma_vap_analytical));
  Modelica.Utilities.Streams.print("  (c_sigma)@vap Span= " + String(c_sigma_vap_Span2000));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("internal Energy");
  // check (du/dT)@liq
  duT_liq_numerical = (sat_Tplus.liq.u-sat_Tminus.liq.u)/(sat_Tplus.liq.T-sat_Tminus.liq.T);
  duT_liq_analytical = duTp_liq+dupT_liq*Medium.saturationPressure_derT(T=T);
  Modelica.Utilities.Streams.print("  (du/dT)@liq  numerical= " + String(duT_liq_numerical));
  Modelica.Utilities.Streams.print("  (du/dT)@liq analytical= " + String(duT_liq_analytical));
  // check (du/dT)@vap
  duT_vap_numerical = (sat_Tplus.vap.u-sat_Tminus.vap.u)/(sat_Tplus.vap.T-sat_Tminus.vap.T);
  duT_vap_analytical = duTp_vap+dupT_vap*Medium.saturationPressure_derT(T=T);
  Modelica.Utilities.Streams.print("  (du/dT)@vap  numerical= " + String(duT_vap_numerical));
  Modelica.Utilities.Streams.print("  (du/dT)@vap analytical= " + String(duT_vap_analytical));

  annotation (experiment(NumberOfIntervals=1));
end Validate_Derivatives_SaturationBoundary;
