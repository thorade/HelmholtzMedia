within HelmholtzMedia.Examples.Validation;
model Derivatives_SaturationBoundary
  "compare analytical derivatives to numerical derivatives"

  package Medium = HelmholtzFluids.Butane;

  // right at T_trip and T_crit, numerical derivatives will fail
  Modelica.Blocks.Sources.Ramp T_ramp(
    duration=10,
    startTime=0,
    height=Tcrit - Tmin - 3,
    offset=Tmin + 1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Medium.Temperature T=T_ramp.y;
  Medium.SaturationProperties sat=Medium.setSat_T(T=T);
  Medium.DerPressureByTemperature dpT=Medium.saturationPressure_derT(T=T);
  Medium.DerTemperatureByPressure dTp=Medium.saturationTemperature_derp(p=sat.psat);

// Density derivatives
  Medium.DerDensityByTemperature ddT_liq_numerical;
  Medium.DerDensityByTemperature ddT_liq_analytical;
  Medium.DerDensityByTemperature ddT_liq_analytical2;
  Medium.DerDensityByTemperature ddT_vap_numerical;
  Medium.DerDensityByTemperature ddT_vap_analytical;
  Medium.DerDensityByPressure ddp_liq_numerical;
  Medium.DerDensityByPressure ddp_liq_analytical1;
  Medium.DerDensityByPressure ddp_liq_analytical2;
  Medium.DerDensityByPressure ddp_vap_numerical;
  Medium.DerDensityByPressure ddp_vap_analytical1;
  Medium.DerDensityByPressure ddp_vap_analytical2;
  Medium.DerVolumeByPressure dvp_liq_analytical1;
  Medium.DerVolumeByPressure dvp_vap_analytical1;
// Enthalpy derivatives
  Medium.DerEnthalpyByTemperature dhT_liq_analytical;
  Medium.DerEnthalpyByTemperature dhT_liq_numerical;
  Medium.DerEnthalpyByTemperature dhT_vap_analytical;
  Medium.DerEnthalpyByTemperature dhT_vap_numerical;
  Medium.DerEnthalpyByPressure dhp_liq_analytical;
  Medium.DerEnthalpyByPressure dhp_liq_numerical;
  Medium.DerEnthalpyByPressure dhp_vap_analytical;
  Medium.DerEnthalpyByPressure dhp_vap_numerical;
// Entropy derivatives
  Medium.DerEntropyByTemperature dsT_liq_analytical;
  Medium.DerEntropyByTemperature dsT_liq_analytical2;
  Medium.DerEntropyByTemperature dsT_liq_numerical;
  Medium.DerEntropyByTemperature dsT_vap_analytical;
  Medium.DerEntropyByTemperature dsT_vap_numerical;
  Medium.DerEntropyByPressure dsp_liq_analytical;
  Medium.DerEntropyByPressure dsp_liq_numerical;
  Medium.DerEntropyByPressure dsp_vap_analytical;
  Medium.DerEntropyByPressure dsp_vap_numerical;
//saturated heat capacity
  Medium.DerEnthalpyByTemperature c_sigma_liq_Span2000;
  Medium.DerEnthalpyByTemperature c_sigma_liq_analytical;
  Medium.DerEnthalpyByTemperature c_sigma_liq_numerical;
  Medium.DerEnthalpyByTemperature c_sigma_vap_Span2000;
  Medium.DerEnthalpyByTemperature c_sigma_vap_analytical;
  Medium.DerEnthalpyByTemperature c_sigma_vap_numerical;
// Energy derivatives
  Medium.DerEnergyByTemperature duT_liq_numerical;
  Medium.DerEnergyByTemperature duT_liq_analytical;
  Medium.DerEnergyByTemperature duT_vap_numerical;
  Medium.DerEnergyByTemperature duT_vap_analytical;

protected
  constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  Medium.EoS.HelmholtzDerivs fl=Medium.EoS.setHelmholtzDerivsSecond(T=T, d=sat.liq.d, phase=1);
  Medium.EoS.HelmholtzDerivs fv=Medium.EoS.setHelmholtzDerivsSecond(T=T, d=sat.vap.d, phase=1);
  Medium.SaturationProperties sat_Tplus = Medium.setSat_T(T=1.0001*T);
  Medium.SaturationProperties sat_Tminus= Medium.setSat_T(T=0.9999*T);
  Medium.SaturationProperties sat_pplus = Medium.setSat_p(p=1.0001*sat.psat);
  Medium.SaturationProperties sat_pminus= Medium.setSat_p(p=0.9999*sat.psat);
// Entropy single phase derivatives
  Medium.DerEntropyByDensity dsdT_liq=fl.R/sat.liq.d*(-(1 + fl.delta*
      fl.rd) + (0 + fl.tau*fl.delta*fl.rtd));
  Medium.DerEntropyByTemperature dsTd_liq=fl.R/T*(-fl.tau^2*(fl.itt
       + fl.rtt));
  Medium.DerEntropyByTemperature dsTp_liq=dsTd_liq - dsdT_liq*
      Medium.pressure_derT_d(state=sat.liq)/Medium.pressure_derd_T(state=sat.liq);
  Medium.DerEntropyByPressure dspT_liq=dsdT_liq/
      Medium.pressure_derd_T(state=sat.liq);
  Medium.DerEntropyByDensity dsdT_vap=fv.R/sat.vap.d*(-(1 + fv.delta*
      fv.rd) + (0 + fv.tau*fv.delta*fv.rtd));
  Medium.DerEntropyByTemperature dsTd_vap=fv.R/T*(-fv.tau^2*(fv.itt
       + fv.rtt));
  Medium.DerEntropyByTemperature dsTp_vap=dsTd_vap - dsdT_vap*
      Medium.pressure_derT_d(state=sat.vap)/Medium.pressure_derd_T(state=sat.vap);
  Medium.DerEntropyByPressure dspT_vap=dsdT_vap/
      Medium.pressure_derd_T(state=sat.vap);
// Internal energy single phase derivatives
  Medium.DerEnergyByDensity dudT_liq=fl.R*T/sat.liq.d*fl.tau*fl.delta
      *fl.rtd;
  Medium.DerEnergyByTemperature duTd_liq=
      Medium.specificHeatCapacityCv(state=sat.liq);
  Medium.DerEnergyByTemperature duTp_liq=duTd_liq - dudT_liq*
      Medium.pressure_derT_d(state=sat.liq)/Medium.pressure_derd_T(state=sat.liq);
  Medium.DerEnergyByPressure dupT_liq=dudT_liq/Medium.pressure_derd_T(
      state=sat.liq);
  Medium.DerEnergyByDensity dudT_vap=fv.R*T/sat.vap.d*fv.tau*fv.delta
      *fv.rtd;
  Medium.DerEnergyByTemperature duTd_vap=
      Medium.specificHeatCapacityCv(state=sat.vap);
  Medium.DerEnergyByTemperature duTp_vap=duTd_vap - dudT_vap*
      Medium.pressure_derT_d(state=sat.vap)/Medium.pressure_derd_T(state=sat.vap);
  Medium.DerEnergyByPressure dupT_vap=dudT_vap/Medium.pressure_derd_T(
      state=sat.vap);

equation
  // Modelica.Utilities.Streams.print(" ");
  // Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|"); // 80 characters
  // Modelica.Utilities.Streams.print(" ");

  // Modelica.Utilities.Streams.print("Density");
  // check (dd/dT)@liq
  ddT_liq_numerical = (sat_Tplus.liq.d - sat_Tminus.liq.d)/(sat_Tplus.liq.T - sat_Tminus.liq.T);
  ddT_liq_analytical = Medium.density_derT_p(state=sat.liq) +Medium.density_derp_T(state=sat.liq)*dpT;
  ddT_liq_analytical2 = (dpT - Medium.EoS.dpTd(fl))/Medium.EoS.dpdT(fl);
  // Modelica.Utilities.Streams.print("  (dd/dT)@liq  numerical= " + String(ddT_liq_numerical));
  // Modelica.Utilities.Streams.print("  (dd/dT)@liq analytical= " + String(ddT_liq_analytical));
  // check (dd/dT)@vap
  ddT_vap_numerical = (sat_Tplus.vap.d - sat_Tminus.vap.d)/(sat_Tplus.vap.T - sat_Tminus.vap.T);
  ddT_vap_analytical = Medium.density_derT_p(state=sat.vap) +Medium.density_derp_T(state=sat.vap)*dpT;
  // Modelica.Utilities.Streams.print("  (dd/dT)@vap  numerical= " + String(ddT_vap_numerical));
  // Modelica.Utilities.Streams.print("  (dd/dT)@vap analytical= " + String(ddT_vap_analytical));
  // check (dd/dp)@liq
  ddp_liq_numerical = (sat_pplus.liq.d - sat_pminus.liq.d)/(sat_pplus.liq.p - sat_pminus.liq.p);
  ddp_liq_analytical1 = Medium.density_derp_T(state=sat.liq) +Medium.density_derT_p(state=sat.liq)*dTp;
  ddp_liq_analytical2 = Medium.dBubbleDensity_dPressure(sat=sat);
  // Modelica.Utilities.Streams.print("  (dd/dp)@liq   numerical= " + String(ddp_liq_numerical));
  // Modelica.Utilities.Streams.print("  (dd/dp)@liq analytical1= " + String(ddp_liq_analytical1));
  // Modelica.Utilities.Streams.print("  (dd/dp)@liq analytical2= " + String(ddp_liq_analytical2));
  // check (dd/dp)@vap
  ddp_vap_numerical = (sat_pplus.vap.d - sat_pminus.vap.d)/(sat_pplus.vap.p - sat_pminus.vap.p);
  ddp_vap_analytical1 = Medium.density_derp_T(state=sat.vap) +Medium.density_derT_p(state=sat.vap)*dTp;
  ddp_vap_analytical2 = Medium.dDewDensity_dPressure(sat=sat);
  // Modelica.Utilities.Streams.print("  (dd/dp)@vap   numerical= " + String(ddp_vap_numerical));
  // Modelica.Utilities.Streams.print("  (dd/dp)@vap analytical1= " + String(ddp_vap_analytical1));
  // Modelica.Utilities.Streams.print("  (dd/dp)@vap analytical2= " + String(ddp_vap_analytical2));
  // check (dv/dp)@liq
  dvp_liq_analytical1 = -1.0/sat.liq.d^2*ddp_liq_analytical1;
  // check (dv/dp)@vap
  dvp_vap_analytical1 = -1.0/sat.vap.d^2*ddp_vap_analytical1;

  // Modelica.Utilities.Streams.print(" ");
  // Modelica.Utilities.Streams.print("Enthalpy");
  // check (dh/dT)@liq
  dhT_liq_numerical = (sat_Tplus.liq.h-sat_Tminus.liq.h)/(sat_Tplus.liq.T-sat_Tminus.liq.T);
  dhT_liq_analytical = Medium.specificHeatCapacityCp(state=sat.liq) +Medium.isothermalThrottlingCoefficient(state=sat.liq)*dpT;
  // Modelica.Utilities.Streams.print("  (dh/dT)@liq  numerical= " + String(dhT_liq_numerical));
  // Modelica.Utilities.Streams.print("  (dh/dT)@liq analytical= " + String(dhT_liq_analytical));
  // check (dh/dT)@vap
  dhT_vap_numerical = (sat_Tplus.vap.h-sat_Tminus.vap.h)/(sat_Tplus.vap.T-sat_Tminus.vap.T);
  dhT_vap_analytical = Medium.specificHeatCapacityCp(state=sat.vap) +Medium.isothermalThrottlingCoefficient(state=sat.vap)*dpT;
  // Modelica.Utilities.Streams.print("  (dh/dT)@vap  numerical= " + String(dhT_vap_numerical));
  // Modelica.Utilities.Streams.print("  (dh/dT)@vap analytical= " + String(dhT_vap_analytical));
  // check (dh/dp)@liq
  dhp_liq_numerical = (sat_pplus.liq.h-sat_pminus.liq.h)/(sat_pplus.liq.p-sat_pminus.liq.p);
  dhp_liq_analytical = Medium.dBubbleEnthalpy_dPressure(sat=sat);
  // Modelica.Utilities.Streams.print("  (dh/dp)@liq  numerical= " + String(dhp_liq_numerical));
  // Modelica.Utilities.Streams.print("  (dh/dp)@liq analytical= " + String(dhp_liq_analytical));
  // check (dh/dp)@vap
  dhp_vap_numerical = (sat_pplus.vap.h-sat_pminus.vap.h)/(sat_pplus.vap.p-sat_pminus.vap.p);
  dhp_vap_analytical = Medium.dDewEnthalpy_dPressure(sat=sat);
  // Modelica.Utilities.Streams.print("  (dh/dp)@vap  numerical= " + String(dhp_vap_numerical));
  // Modelica.Utilities.Streams.print("  (dh/dp)@vap analytical= " + String(dhp_vap_analytical));

  // Modelica.Utilities.Streams.print(" ");
  // Modelica.Utilities.Streams.print("Entropy");
  // check (ds/dT)@liq
  dsT_liq_numerical = (sat_Tplus.liq.s-sat_Tminus.liq.s)/(sat_Tplus.liq.T-sat_Tminus.liq.T);
  dsT_liq_analytical = dsTp_liq+dspT_liq*dpT;
  dsT_liq_analytical2 = Medium.EoS.dsTd(fl) + Medium.EoS.dsdT(fl)*(dpT - Medium.EoS.dpTd(fl))/Medium.EoS.dpdT(fl);
  // Modelica.Utilities.Streams.print("  (ds/dT)@liq  numerical= " + String(dsT_liq_numerical));
  // Modelica.Utilities.Streams.print("  (ds/dT)@liq analytical= " + String(dsT_liq_analytical));
  // check (ds/dT)@vap
  dsT_vap_numerical = (sat_Tplus.vap.s-sat_Tminus.vap.s)/(sat_Tplus.vap.T-sat_Tminus.vap.T);
  dsT_vap_analytical = dsTp_vap+dspT_vap*dpT;
  // Modelica.Utilities.Streams.print("  (ds/dT)@vap  numerical= " + String(dsT_vap_numerical));
  // Modelica.Utilities.Streams.print("  (ds/dT)@vap analytical= " + String(dsT_vap_analytical));
  // check (ds/dp)@liq
  dsp_liq_numerical = (sat_pplus.liq.s-sat_pminus.liq.s)/(sat_pplus.liq.p-sat_pminus.liq.p);
  dsp_liq_analytical = dspT_liq+dsTp_liq*dTp;
  // Modelica.Utilities.Streams.print("  (ds/dp)@liq  numerical= " + String(dsp_liq_numerical));
  // Modelica.Utilities.Streams.print("  (ds/dp)@liq analytical= " + String(dsp_liq_analytical));
  // check (ds/dp)@vap
  dsp_vap_numerical = (sat_pplus.vap.s-sat_pminus.vap.s)/(sat_pplus.vap.p-sat_pminus.vap.p);
  dsp_vap_analytical = dspT_vap+dsTp_vap*dTp;
  // Modelica.Utilities.Streams.print("  (ds/dp)@vap  numerical= " + String(dsp_vap_numerical));
  // Modelica.Utilities.Streams.print("  (ds/dp)@vap analytical= " + String(dsp_vap_analytical));

  // Modelica.Utilities.Streams.print(" ");
  // Modelica.Utilities.Streams.print("saturated heat capacity");
  // check (c_sigma)@liq
  c_sigma_liq_numerical = T*dsT_liq_numerical;
  c_sigma_liq_analytical = T*dsT_liq_analytical;
  c_sigma_liq_Span2000 = Medium.specificHeatCapacityCv(state=sat.liq) - T*Medium.pressure_derT_d(state=sat.liq)/(sat.liq.d^2)*ddT_liq_analytical;
  // Modelica.Utilities.Streams.print("  (c_sigma)@liq numerical= " + String(c_sigma_liq_numerical));
  // Modelica.Utilities.Streams.print("  (c_sigma)@liq analytical= " + String(c_sigma_liq_analytical));
  // Modelica.Utilities.Streams.print("  (c_sigma)@liq Span= " + String(c_sigma_liq_Span2000));
  // check (c_sigma)@vap
  c_sigma_vap_numerical = T*dsT_vap_numerical;
  c_sigma_vap_analytical = T*dsT_vap_analytical;
  c_sigma_vap_Span2000 = Medium.specificHeatCapacityCv(state=sat.vap) - T*Medium.pressure_derT_d(state=sat.vap)/(sat.vap.d^2)*ddT_vap_analytical;
  // Modelica.Utilities.Streams.print("  (c_sigma)@vap numerical= " + String(c_sigma_vap_numerical));
  // Modelica.Utilities.Streams.print("  (c_sigma)@vap analytical= " + String(c_sigma_vap_analytical));
  // Modelica.Utilities.Streams.print("  (c_sigma)@vap Span= " + String(c_sigma_vap_Span2000));

  // Modelica.Utilities.Streams.print(" ");
  // Modelica.Utilities.Streams.print("internal Energy");
  // check (du/dT)@liq
  duT_liq_numerical = (sat_Tplus.liq.u-sat_Tminus.liq.u)/(sat_Tplus.liq.T-sat_Tminus.liq.T);
  duT_liq_analytical = duTp_liq+dupT_liq*dpT;
  // Modelica.Utilities.Streams.print("  (du/dT)@liq  numerical= " + String(duT_liq_numerical));
  // Modelica.Utilities.Streams.print("  (du/dT)@liq analytical= " + String(duT_liq_analytical));
  // check (du/dT)@vap
  duT_vap_numerical = (sat_Tplus.vap.u-sat_Tminus.vap.u)/(sat_Tplus.vap.T-sat_Tminus.vap.T);
  duT_vap_analytical = duTp_vap+dupT_vap*dpT;
  // Modelica.Utilities.Streams.print("  (du/dT)@vap  numerical= " + String(duT_vap_numerical));
  // Modelica.Utilities.Streams.print("  (du/dT)@vap analytical= " + String(duT_vap_analytical));

  annotation (experiment(StopTime=10));
end Derivatives_SaturationBoundary;
