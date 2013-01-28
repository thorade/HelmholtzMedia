within HelmholtzMedia.Examples.Validation;
model Derivatives_Helmholtz_Multiplied
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given by RefProp
  // go to Options, Preferences, check "Show options used for analyzing EoS"

  package medium = HelmholtzMedia.HelmholtzFluids.Propane;
  parameter medium.Density d=1e-5;
  parameter medium.Temperature T=500;

protected
  String fileName = "HelmholtzDerivs_multiplied.csv";
  // While csv originally stood for comma-seperated-values, MS Excel uses semicolons to seperate the values
  String Separator = ";";

  constant medium.MolarMass MM = medium.fluidConstants[1].molarMass;
  constant medium.SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant medium.Density d_crit=MM/medium.fluidConstants[1].criticalMolarVolume;
  constant medium.Temperature T_crit=medium.fluidConstants[1].criticalTemperature;
  constant medium.Temperature T_trip=medium.fluidConstants[1].triplePointTemperature;

  medium.SaturationProperties sat_trip = medium.setSat_T(T=T_trip);
  medium.SaturationProperties sat_IIR = medium.setSat_T(T=273.15); // 0°C
  medium.SaturationProperties sat_ASHRAE = medium.setSat_T(T=233.15); // -40°C
  medium.SaturationProperties sat_NBP = medium.setSat_p(p=101325); // 1.01325 bar = 1atm

  medium.EoS.HelmholtzDerivs f_crit = medium.EoS.setHelmholtzDerivsThird(T=T_crit, d=d_crit, phase=1);
  medium.EoS.HelmholtzDerivs f_tl = medium.EoS.setHelmholtzDerivsThird(T=sat_trip.liq.T, d=sat_trip.liq.d, phase=1);
  medium.EoS.HelmholtzDerivs f_tv = medium.EoS.setHelmholtzDerivsThird(T=sat_trip.vap.T, d=sat_trip.vap.d, phase=1);
  medium.EoS.HelmholtzDerivs f_IIR = medium.EoS.setHelmholtzDerivsThird(T=sat_IIR.liq.T, d=sat_IIR.liq.d, phase=1);
  medium.EoS.HelmholtzDerivs f_ASHRAE = medium.EoS.setHelmholtzDerivsThird(T=sat_ASHRAE.liq.T, d=sat_ASHRAE.liq.d, phase=1);
  medium.EoS.HelmholtzDerivs f_NBP = medium.EoS.setHelmholtzDerivsThird(T=sat_NBP.liq.T, d=sat_NBP.liq.d, phase=1);
  medium.EoS.HelmholtzDerivs f = medium.EoS.setHelmholtzDerivsThird(T=T, d=d, phase=1);
  medium.EoS.HelmholtzDerivs f_num(T=T, d=d);

  Real delta(unit="1")=d/d_crit "reduced density";
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";
  Real eps= 1e-6;

algorithm
  // numerical derivative, for comparison, last two line of csv file should be identical
  f_num.i    := medium.EoS.f_i(tau=tau, delta=delta);
  f_num.it   := (medium.EoS.f_i(tau=tau+eps, delta=delta)-medium.EoS.f_i(tau=tau-eps, delta=delta))/(2*eps);
  f_num.itt  := (medium.EoS.f_it(tau=tau+eps, delta=delta)-medium.EoS.f_it(tau=tau-eps, delta=delta))/(2*eps);
  f_num.ittt  := (medium.EoS.f_itt(tau=tau+eps, delta=delta)-medium.EoS.f_itt(tau=tau-eps, delta=delta))/(2*eps);

  f_num.r    := medium.EoS.f_r(tau=tau, delta=delta);
  f_num.rt   := (medium.EoS.f_r(tau=tau+eps, delta=delta)-medium.EoS.f_r(tau=tau-eps, delta=delta))/(2*eps);
  f_num.rtt  := (medium.EoS.f_rt(tau=tau+eps, delta=delta)-medium.EoS.f_rt(tau=tau-eps, delta=delta))/(2*eps);
  f_num.rttt  := (medium.EoS.f_rtt(tau=tau+eps, delta=delta)-medium.EoS.f_rtt(tau=tau-eps, delta=delta))/(2*eps);
  f_num.rtd  := (medium.EoS.f_rt(tau=tau, delta=delta+eps)-medium.EoS.f_rt(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rttd := (medium.EoS.f_rtt(tau=tau, delta=delta+eps)-medium.EoS.f_rtt(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rtdd :=(medium.EoS.f_rtd(tau=tau, delta=delta+eps)-medium.EoS.f_rtd(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rd   := (medium.EoS.f_r(tau=tau, delta=delta+eps)-medium.EoS.f_r(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rdd  := (medium.EoS.f_rd(tau=tau, delta=delta+eps)-medium.EoS.f_rd(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rddd := (medium.EoS.f_rdd(tau=tau, delta=delta+eps)-medium.EoS.f_rdd(tau=tau, delta=delta-eps))/(2*eps);

  // remove old file
  Modelica.Utilities.Files.remove(fileName);

  // print headers
  Modelica.Utilities.Streams.print("T" +Separator
                                 + "d" +Separator
                                 + "tau" +Separator
                                 + "delta" +Separator
                                 + "alpha_i" +Separator
                                 + "tau*alpha_it" +Separator
                                 + "tau*tau*alpha_itt" +Separator
                                 + "tau*tau*tau*alpha_ittt" +Separator
                                 + "alpha_r" +Separator
                                 + "tau*alpha_rt" +Separator
                                 + "tau*tau*alpha_rtt" +Separator
                                 + "tau*tau*tau*alpha_rttt" +Separator
                                 + "delta*alpha_rd" +Separator
                                 + "delta*delta*alpha_rdd" +Separator
                                 + "delta*delta*delta*alpha_rddd" +Separator
                                 + "tau*delta*alpha_rtd" +Separator
                                 + "tau*delta*delta*alpha_rtdd" +Separator
                                 + "tau*tau*delta*alpha_rttd" +Separator,
                                   fileName);

  // print the actual values
  Modelica.Utilities.Streams.print(String(f_crit.T) + Separator
                                 + String(f_crit.d) + Separator
                                 + String(f_crit.tau) + Separator
                                 + String(f_crit.delta) + Separator
                                 + String(f_crit.i) + Separator
                                 + String(f_crit.it*f_crit.tau)+Separator
                                 + String(f_crit.itt*f_crit.tau*f_crit.tau)+Separator
                                 + String(f_crit.ittt*f_crit.tau*f_crit.tau*f_crit.tau)+Separator
                                 + String(f_crit.r)+Separator
                                 + String(f_crit.rt*f_crit.tau)+Separator
                                 + String(f_crit.rtt*f_crit.tau*f_crit.tau)+Separator
                                 + String(f_crit.rttt*f_crit.tau*f_crit.tau*f_crit.tau)+Separator
                                 + String(f_crit.rd*f_crit.delta)+Separator
                                 + String(f_crit.rdd*f_crit.delta*f_crit.delta)+Separator
                                 + String(f_crit.rddd*f_crit.delta*f_crit.delta*f_crit.delta)+Separator
                                 + String(f_crit.rtd*f_crit.tau*f_crit.delta)+Separator
                                 + String(f_crit.rtdd*f_crit.tau*f_crit.delta*f_crit.delta)+Separator
                                 + String(f_crit.rttd*f_crit.tau*f_crit.tau*f_crit.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_tl.T) + Separator
                                 + String(f_tl.d) + Separator
                                 + String(f_tl.tau) + Separator
                                 + String(f_tl.delta) + Separator
                                 + String(f_tl.i) + Separator
                                 + String(f_tl.it*f_tl.tau)+Separator
                                 + String(f_tl.itt*f_tl.tau*f_tl.tau)+Separator
                                 + String(f_tl.ittt*f_tl.tau*f_tl.tau*f_tl.tau)+Separator
                                 + String(f_tl.r)+Separator
                                 + String(f_tl.rt*f_tl.tau)+Separator
                                 + String(f_tl.rtt*f_tl.tau*f_tl.tau)+Separator
                                 + String(f_tl.rttt*f_tl.tau*f_tl.tau*f_tl.tau)+Separator
                                 + String(f_tl.rd*f_tl.delta)+Separator
                                 + String(f_tl.rdd*f_tl.delta*f_tl.delta)+Separator
                                 + String(f_tl.rddd*f_tl.delta*f_tl.delta*f_tl.delta)+Separator
                                 + String(f_tl.rtd*f_tl.tau*f_tl.delta)+Separator
                                 + String(f_tl.rtdd*f_tl.tau*f_tl.delta*f_tl.delta)+Separator
                                 + String(f_tl.rttd*f_tl.tau*f_tl.tau*f_tl.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_tv.T) + Separator
                                 + String(f_tv.d) + Separator
                                 + String(f_tv.tau) + Separator
                                 + String(f_tv.delta) + Separator
                                 + String(f_tv.i) + Separator
                                 + String(f_tv.it*f_tv.tau)+Separator
                                 + String(f_tv.itt*f_tv.tau*f_tv.tau)+Separator
                                 + String(f_tv.ittt*f_tv.tau*f_tv.tau*f_tv.tau)+Separator
                                 + String(f_tv.r)+Separator
                                 + String(f_tv.rt*f_tv.tau)+Separator
                                 + String(f_tv.rtt*f_tv.tau*f_tv.tau)+Separator
                                 + String(f_tv.rttt*f_tv.tau*f_tv.tau*f_tv.tau)+Separator
                                 + String(f_tv.rd*f_tv.delta)+Separator
                                 + String(f_tv.rdd*f_tv.delta*f_tv.delta)+Separator
                                 + String(f_tv.rddd*f_tv.delta*f_tv.delta*f_tv.delta)+Separator
                                 + String(f_tv.rtd*f_tv.tau*f_tv.delta)+Separator
                                 + String(f_tv.rtdd*f_tv.tau*f_tv.delta*f_tv.delta)+Separator
                                 + String(f_tv.rttd*f_tv.tau*f_tv.tau*f_tv.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_IIR.T) + Separator
                                 + String(f_IIR.d) + Separator
                                 + String(f_IIR.tau) + Separator
                                 + String(f_IIR.delta) + Separator
                                 + String(f_IIR.i) + Separator
                                 + String(f_IIR.it*f_IIR.tau)+Separator
                                 + String(f_IIR.itt*f_IIR.tau*f_IIR.tau)+Separator
                                 + String(f_IIR.ittt*f_IIR.tau*f_IIR.tau*f_IIR.tau)+Separator
                                 + String(f_IIR.r)+Separator
                                 + String(f_IIR.rt*f_IIR.tau)+Separator
                                 + String(f_IIR.rtt*f_IIR.tau*f_IIR.tau)+Separator
                                 + String(f_IIR.rttt*f_IIR.tau*f_IIR.tau*f_IIR.tau)+Separator
                                 + String(f_IIR.rd*f_IIR.delta)+Separator
                                 + String(f_IIR.rdd*f_IIR.delta*f_IIR.delta)+Separator
                                 + String(f_IIR.rddd*f_IIR.delta*f_IIR.delta*f_IIR.delta)+Separator
                                 + String(f_IIR.rtd*f_IIR.tau*f_IIR.delta)+Separator
                                 + String(f_IIR.rtdd*f_IIR.tau*f_IIR.delta*f_IIR.delta)+Separator
                                 + String(f_IIR.rttd*f_IIR.tau*f_IIR.tau*f_IIR.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_ASHRAE.T) + Separator
                                 + String(f_ASHRAE.d) + Separator
                                 + String(f_ASHRAE.tau) + Separator
                                 + String(f_ASHRAE.delta) + Separator
                                 + String(f_ASHRAE.i) + Separator
                                 + String(f_ASHRAE.it*f_ASHRAE.tau)+Separator
                                 + String(f_ASHRAE.itt*f_ASHRAE.tau*f_ASHRAE.tau)+Separator
                                 + String(f_ASHRAE.ittt*f_ASHRAE.tau*f_ASHRAE.tau*f_ASHRAE.tau)+Separator
                                 + String(f_ASHRAE.r)+Separator
                                 + String(f_ASHRAE.rt*f_ASHRAE.tau)+Separator
                                 + String(f_ASHRAE.rtt*f_ASHRAE.tau*f_ASHRAE.tau)+Separator
                                 + String(f_ASHRAE.rttt*f_ASHRAE.tau*f_ASHRAE.tau*f_ASHRAE.tau)+Separator
                                 + String(f_ASHRAE.rd*f_ASHRAE.delta)+Separator
                                 + String(f_ASHRAE.rdd*f_ASHRAE.delta*f_ASHRAE.delta)+Separator
                                 + String(f_ASHRAE.rddd*f_ASHRAE.delta*f_ASHRAE.delta*f_ASHRAE.delta)+Separator
                                 + String(f_ASHRAE.rtd*f_ASHRAE.tau*f_ASHRAE.delta)+Separator
                                 + String(f_ASHRAE.rtdd*f_ASHRAE.tau*f_ASHRAE.delta*f_ASHRAE.delta)+Separator
                                 + String(f_ASHRAE.rttd*f_ASHRAE.tau*f_ASHRAE.tau*f_ASHRAE.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_NBP.T) + Separator
                                 + String(f_NBP.d) + Separator
                                 + String(f_NBP.tau) + Separator
                                 + String(f_NBP.delta) + Separator
                                 + String(f_NBP.i) + Separator
                                 + String(f_NBP.it*f_NBP.tau)+Separator
                                 + String(f_NBP.itt*f_NBP.tau*f_NBP.tau)+Separator
                                 + String(f_NBP.ittt*f_NBP.tau*f_NBP.tau*f_NBP.tau)+Separator
                                 + String(f_NBP.r)+Separator
                                 + String(f_NBP.rt*f_NBP.tau)+Separator
                                 + String(f_NBP.rtt*f_NBP.tau*f_NBP.tau)+Separator
                                 + String(f_NBP.rttt*f_NBP.tau*f_NBP.tau*f_NBP.tau)+Separator
                                 + String(f_NBP.rd*f_NBP.delta)+Separator
                                 + String(f_NBP.rdd*f_NBP.delta*f_NBP.delta)+Separator
                                 + String(f_NBP.rddd*f_NBP.delta*f_NBP.delta*f_NBP.delta)+Separator
                                 + String(f_NBP.rtd*f_NBP.tau*f_NBP.delta)+Separator
                                 + String(f_NBP.rtdd*f_NBP.tau*f_NBP.delta*f_NBP.delta)+Separator
                                 + String(f_NBP.rttd*f_NBP.tau*f_NBP.tau*f_NBP.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f.T) + Separator
                                 + String(f.d) + Separator
                                 + String(f.tau) + Separator
                                 + String(f.delta) + Separator
                                 + String(f.i) + Separator
                                 + String(f.it*f.tau)+Separator
                                 + String(f.itt*f.tau*f.tau)+Separator
                                 + String(f.ittt*f.tau*f.tau*f.tau)+Separator
                                 + String(f.r)+Separator
                                 + String(f.rt*f.tau)+Separator
                                 + String(f.rtt*f.tau*f.tau)+Separator
                                 + String(f.rttt*f.tau*f.tau*f.tau)+Separator
                                 + String(f.rd*f.delta)+Separator
                                 + String(f.rdd*f.delta*f.delta)+Separator
                                 + String(f.rddd*f.delta*f.delta*f.delta)+Separator
                                 + String(f.rtd*f.tau*f.delta)+Separator
                                 + String(f.rtdd*f.tau*f.delta*f.delta)+Separator
                                 + String(f.rttd*f.tau*f.tau*f.delta)+Separator,
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_num.T) + Separator
                                 + String(f_num.d) + Separator
                                 + String(f_num.tau) + Separator
                                 + String(f_num.delta) + Separator
                                 + String(f_num.i) + Separator
                                 + String(f_num.it*f_num.tau)+Separator
                                 + String(f_num.itt*f_num.tau*f_num.tau)+Separator
                                 + String(f_num.ittt*f_num.tau*f_num.tau*f_num.tau)+Separator
                                 + String(f_num.r)+Separator
                                 + String(f_num.rt*f_num.tau)+Separator
                                 + String(f_num.rtt*f_num.tau*f_num.tau)+Separator
                                 + String(f_num.rttt*f_num.tau*f_num.tau*f_num.tau)+Separator
                                 + String(f_num.rd*f_num.delta)+Separator
                                 + String(f_num.rdd*f_num.delta*f_num.delta)+Separator
                                 + String(f_num.rddd*f_num.delta*f_num.delta*f_num.delta)+Separator
                                 + String(f_num.rtd*f_num.tau*f_num.delta)+Separator
                                 + String(f_num.rtdd*f_num.tau*f_num.delta*f_num.delta)+Separator
                                 + String(f_num.rttd*f_num.tau*f_num.tau*f_num.delta)+Separator,
                                   fileName);

annotation (experiment(NumberOfIntervals=1));
end Derivatives_Helmholtz_Multiplied;
