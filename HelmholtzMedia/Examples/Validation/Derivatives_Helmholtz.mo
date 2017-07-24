within HelmholtzMedia.Examples.Validation;
model Derivatives_Helmholtz
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide;
  parameter Medium.Density d=533;
  parameter Medium.Temperature T=304;
  Medium.ThermodynamicState state = Medium.setState_dTX(d=d, T=T, phase=1);
  Medium.EoS.HelmholtzDerivs f = Medium.EoS.setHelmholtzDerivsThird(T=T, d=d, phase=1);

protected
  String fileName = "HelmholtzDerivs.csv";
  // While csv originally stood for comma-seperated-values, MS Excel uses semicolons to seperate the values
  String Separator = ";";

  constant Medium.MolarMass MM = Medium.fluidConstants[1].molarMass;
  constant Medium.SpecificHeatCapacity R=Medium.fluidConstants[1].gasConstant/MM
    "specific gas constant";
  constant Medium.Density d_crit=MM/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Temperature T_trip=Medium.fluidConstants[1].triplePointTemperature;

  Medium.SaturationProperties sat_trip = Medium.setSat_T(T=T_trip);
  Medium.SaturationProperties sat_IIR = Medium.setSat_T(T=273.15); // 0°C
  Medium.SaturationProperties sat_ASHRAE = Medium.setSat_T(T=233.15); // -40°C
  Medium.SaturationProperties sat_NBP = Medium.setSat_p(p=101325); // 1.01325 bar = 1atm

  Medium.EoS.HelmholtzDerivs f_crit = Medium.EoS.setHelmholtzDerivsThird(T=T_crit, d=d_crit, phase=1);
  Medium.EoS.HelmholtzDerivs f_tl = Medium.EoS.setHelmholtzDerivsThird(T=sat_trip.liq.T, d=sat_trip.liq.d, phase=1);
  Medium.EoS.HelmholtzDerivs f_tv = Medium.EoS.setHelmholtzDerivsThird(T=sat_trip.vap.T, d=sat_trip.vap.d, phase=1);
  Medium.EoS.HelmholtzDerivs f_IIR = Medium.EoS.setHelmholtzDerivsThird(T=sat_IIR.liq.T, d=sat_IIR.liq.d, phase=1);
  Medium.EoS.HelmholtzDerivs f_ASHRAE = Medium.EoS.setHelmholtzDerivsThird(T=sat_ASHRAE.liq.T, d=sat_ASHRAE.liq.d, phase=1);
  Medium.EoS.HelmholtzDerivs f_NBP = Medium.EoS.setHelmholtzDerivsThird(T=sat_NBP.liq.T, d=sat_NBP.liq.d, phase=1);
  Medium.EoS.HelmholtzDerivs f_num(T=T, d=d);

  Real delta(unit="1")=d/d_crit "reduced density";
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";
  Real eps= 1e-6;

algorithm
  // numerical derivative, for comparison, last two line of csv file should be identical
  f_num.i    := Medium.EoS.f_i(tau=tau, delta=delta);
  f_num.it   := (Medium.EoS.f_i(tau=tau+eps, delta=delta)-Medium.EoS.f_i(tau=tau-eps, delta=delta))/(2*eps);
  f_num.itt  := (Medium.EoS.f_it(tau=tau+eps, delta=delta)-Medium.EoS.f_it(tau=tau-eps, delta=delta))/(2*eps);
  f_num.ittt  := (Medium.EoS.f_itt(tau=tau+eps, delta=delta)-Medium.EoS.f_itt(tau=tau-eps, delta=delta))/(2*eps);

  f_num.r    := Medium.EoS.f_r(tau=tau, delta=delta);
  f_num.rt   := (Medium.EoS.f_r(tau=tau+eps, delta=delta)-Medium.EoS.f_r(tau=tau-eps, delta=delta))/(2*eps);
  f_num.rtt  := (Medium.EoS.f_rt(tau=tau+eps, delta=delta)-Medium.EoS.f_rt(tau=tau-eps, delta=delta))/(2*eps);
  f_num.rttt  := (Medium.EoS.f_rtt(tau=tau+eps, delta=delta)-Medium.EoS.f_rtt(tau=tau-eps, delta=delta))/(2*eps);
  f_num.rtd  := (Medium.EoS.f_rt(tau=tau, delta=delta+eps)-Medium.EoS.f_rt(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rttd := (Medium.EoS.f_rtt(tau=tau, delta=delta+eps)-Medium.EoS.f_rtt(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rtdd :=(Medium.EoS.f_rtd(tau=tau, delta=delta+eps)-Medium.EoS.f_rtd(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rd   := (Medium.EoS.f_r(tau=tau, delta=delta+eps)-Medium.EoS.f_r(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rdd  := (Medium.EoS.f_rd(tau=tau, delta=delta+eps)-Medium.EoS.f_rd(tau=tau, delta=delta-eps))/(2*eps);
  f_num.rddd := (Medium.EoS.f_rdd(tau=tau, delta=delta+eps)-Medium.EoS.f_rdd(tau=tau, delta=delta-eps))/(2*eps);

  when initial() then
  if not Modelica.Utilities.Files.exist(fileName) then
    // print headers
    Modelica.Utilities.Streams.print("T" +Separator
                                   + "d" +Separator
                                   + "tau" +Separator
                                   + "delta" +Separator
                                   + "alpha_i" +Separator
                                   + "alpha_it" +Separator
                                   + "alpha_itt" +Separator
                                   + "alpha_ittt" +Separator
                                   + "alpha_r" +Separator
                                   + "alpha_rd" +Separator
                                   + "alpha_rdd" +Separator
                                   + "alpha_rt" +Separator
                                   + "alpha_rtt" +Separator
                                   + "alpha_rtd" +Separator
                                   + "alpha_rddd" +Separator
                                   + "alpha_rtdd" +Separator
                                   + "alpha_rttd" +Separator
                                   + "alpha_rttt" +Separator,
                                     fileName);
    // print fixed values
     Modelica.Utilities.Streams.print(String(f_crit.T) + Separator
                                    + String(f_crit.d) + Separator
                                    + String(f_crit.tau) + Separator
                                    + String(f_crit.delta) + Separator
                                    + String(f_crit.i) + Separator
                                    + String(f_crit.it)+Separator
                                    + String(f_crit.itt)+Separator
                                    + String(f_crit.ittt)+Separator
                                    + String(f_crit.r)+Separator
                                    + String(f_crit.rd)+Separator
                                    + String(f_crit.rdd)+Separator
                                    + String(f_crit.rt)+Separator
                                    + String(f_crit.rtt)+Separator
                                    + String(f_crit.rtd)+Separator
                                    + String(f_crit.rddd)+Separator
                                    + String(f_crit.rtdd)+Separator
                                    + String(f_crit.rttd)+Separator
                                    + String(f_crit.rttt)+Separator,
                                      fileName);
    Modelica.Utilities.Streams.print(String(f_tl.T) + Separator
                                   + String(f_tl.d) + Separator
                                   + String(f_tl.tau) + Separator
                                   + String(f_tl.delta) + Separator
                                   + String(f_tl.i) + Separator
                                   + String(f_tl.it)+Separator
                                   + String(f_tl.itt)+Separator
                                   + String(f_tl.ittt)+Separator
                                   + String(f_tl.r)+Separator
                                   + String(f_tl.rd)+Separator
                                   + String(f_tl.rdd)+Separator
                                   + String(f_tl.rt)+Separator
                                   + String(f_tl.rtt)+Separator
                                   + String(f_tl.rtd)+Separator
                                   + String(f_tl.rddd)+Separator
                                   + String(f_tl.rtdd)+Separator
                                   + String(f_tl.rttd)+Separator
                                   + String(f_tl.rttt)+Separator,
                                     fileName);
    Modelica.Utilities.Streams.print(String(f_tv.T) + Separator
                                   + String(f_tv.d) + Separator
                                   + String(f_tv.tau) + Separator
                                   + String(f_tv.delta) + Separator
                                   + String(f_tv.i) + Separator
                                   + String(f_tv.it)+Separator
                                   + String(f_tv.itt)+Separator
                                   + String(f_tv.ittt)+Separator
                                   + String(f_tv.r)+Separator
                                   + String(f_tv.rd)+Separator
                                   + String(f_tv.rdd)+Separator
                                   + String(f_tv.rt)+Separator
                                   + String(f_tv.rtt)+Separator
                                   + String(f_tv.rtd)+Separator
                                   + String(f_tv.rddd)+Separator
                                   + String(f_tv.rtdd)+Separator
                                   + String(f_tv.rttd)+Separator
                                   + String(f_tv.rttt)+Separator,
                                     fileName);
    Modelica.Utilities.Streams.print(String(f_IIR.T) + Separator
                                   + String(f_IIR.d) + Separator
                                   + String(f_IIR.tau) + Separator
                                   + String(f_IIR.delta) + Separator
                                   + String(f_IIR.i) + Separator
                                   + String(f_IIR.it)+Separator
                                   + String(f_IIR.itt)+Separator
                                   + String(f_IIR.ittt)+Separator
                                   + String(f_IIR.r)+Separator
                                   + String(f_IIR.rd)+Separator
                                   + String(f_IIR.rdd)+Separator
                                   + String(f_IIR.rt)+Separator
                                   + String(f_IIR.rtt)+Separator
                                   + String(f_IIR.rtd)+Separator
                                   + String(f_IIR.rddd)+Separator
                                   + String(f_IIR.rtdd)+Separator
                                   + String(f_IIR.rttd)+Separator
                                   + String(f_IIR.rttt)+Separator,
                                     fileName);
    Modelica.Utilities.Streams.print(String(f_ASHRAE.T) + Separator
                                   + String(f_ASHRAE.d) + Separator
                                   + String(f_ASHRAE.tau) + Separator
                                   + String(f_ASHRAE.delta) + Separator
                                   + String(f_ASHRAE.i) + Separator
                                   + String(f_ASHRAE.it)+Separator
                                   + String(f_ASHRAE.itt)+Separator
                                   + String(f_ASHRAE.ittt)+Separator
                                   + String(f_ASHRAE.r)+Separator
                                   + String(f_ASHRAE.rd)+Separator
                                   + String(f_ASHRAE.rdd)+Separator
                                   + String(f_ASHRAE.rt)+Separator
                                   + String(f_ASHRAE.rtt)+Separator
                                   + String(f_ASHRAE.rtd)+Separator
                                   + String(f_ASHRAE.rddd)+Separator
                                   + String(f_ASHRAE.rtdd)+Separator
                                   + String(f_ASHRAE.rttd)+Separator
                                   + String(f_ASHRAE.rttt)+Separator,
                                     fileName);
    Modelica.Utilities.Streams.print(String(f_NBP.T) + Separator
                                   + String(f_NBP.d) + Separator
                                   + String(f_NBP.tau) + Separator
                                   + String(f_NBP.delta) + Separator
                                   + String(f_NBP.i) + Separator
                                   + String(f_NBP.it)+Separator
                                   + String(f_NBP.itt)+Separator
                                   + String(f_NBP.ittt)+Separator
                                   + String(f_NBP.r)+Separator
                                   + String(f_NBP.rd)+Separator
                                   + String(f_NBP.rdd)+Separator
                                   + String(f_NBP.rt)+Separator
                                   + String(f_NBP.rtt)+Separator
                                   + String(f_NBP.rtd)+Separator
                                   + String(f_NBP.rddd)+Separator
                                   + String(f_NBP.rtdd)+Separator
                                   + String(f_NBP.rttd)+Separator
                                   + String(f_NBP.rttt)+Separator,
                                     fileName);
  end if;
  end when;

  when terminal() then
    // print values
    Modelica.Utilities.Streams.print(String(f.T) + Separator
                                   + String(f.d) + Separator
                                   + String(f.tau) + Separator
                                   + String(f.delta) + Separator
                                   + String(f.i) + Separator
                                   + String(f.it)+Separator
                                   + String(f.itt)+Separator
                                   + String(f.ittt)+Separator
                                   + String(f.r)+Separator
                                   + String(f.rd)+Separator
                                   + String(f.rdd)+Separator
                                   + String(f.rt)+Separator
                                   + String(f.rtt)+Separator
                                   + String(f.rtd)+Separator
                                   + String(f.rddd)+Separator
                                   + String(f.rtdd)+Separator
                                   + String(f.rttd)+Separator
                                   + String(f.rttt)+Separator,
                                     fileName);
    Modelica.Utilities.Streams.print(String(f_num.T) + Separator
                                   + String(f_num.d) + Separator
                                   + String(f_num.tau) + Separator
                                   + String(f_num.delta) + Separator
                                   + String(f_num.i) + Separator
                                   + String(f_num.it)+Separator
                                   + String(f_num.itt)+Separator
                                   + String(f_num.ittt)+Separator
                                   + String(f_num.r)+Separator
                                   + String(f_num.rd)+Separator
                                   + String(f_num.rdd)+Separator
                                   + String(f_num.rt)+Separator
                                   + String(f_num.rtt)+Separator
                                   + String(f_num.rtd)+Separator
                                   + String(f_num.rddd)+Separator
                                   + String(f_num.rtdd)+Separator
                                   + String(f_num.rttd)+Separator
                                   + String(f_num.rttt)+Separator,
                                     fileName);
  end when;

  annotation (experiment(Interval=1));
end Derivatives_Helmholtz;
