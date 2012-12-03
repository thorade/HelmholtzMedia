within HelmholtzMedia.Examples.Validation;
model printHelmholtzDerivatives
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  package medium = HelmholtzMedia.HelmholtzFluids.Butane;
  parameter medium.Density d=838.025;
  parameter medium.Temperature T=500;

protected
  String fileName = "Helmholtz_energy_derivatives.csv";

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

algorithm
  // While csv originally stood for comma-seperated-values, MS Excel uses semicolons to seperate the values
  // remove old file
  Modelica.Utilities.Files.remove(fileName);
  // fluid name or CAS number
  Modelica.Utilities.Streams.print("CAS registry number;" + medium.fluidConstants[1].casRegistryNumber,fileName);
  // print headers
  Modelica.Utilities.Streams.print(";;;;ideal;;;;residual;", fileName);
  Modelica.Utilities.Streams.print("T;d;tau;delta;alpha_i;alpha_it;alpha_itt;alpha_ittt;alpha_r;alpha_rd;alpha_rdd;alpha_rt;alpha_rtt;alpha_rtd;alpha_rddd;alpha_rtdd;alpha_rttd;alpha_rttt", fileName);
  // the actual values
  Modelica.Utilities.Streams.print(String(f_crit.T) + ";"
                                 + String(f_crit.d) + ";"
                                 + String(f_crit.tau) + ";"
                                 + String(f_crit.delta) + ";"
                                 + String(f_crit.i) + ";"
                                 + String(f_crit.it)+";"
                                 + String(f_crit.itt)+";"
                                 + String(f_crit.ittt)+";"
                                 + String(f_crit.r)+";"
                                 + String(f_crit.rd)+";"
                                 + String(f_crit.rdd)+";"
                                 + String(f_crit.rt)+";"
                                 + String(f_crit.rtt)+";"
                                 + String(f_crit.rtd)+";"
                                 + String(f_crit.rddd)+";"
                                 + String(f_crit.rtdd)+";"
                                 + String(f_crit.rttd)+";"
                                 + String(f_crit.rttt)+";",
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_tl.T) + ";"
                                 + String(f_tl.d) + ";"
                                 + String(f_tl.tau) + ";"
                                 + String(f_tl.delta) + ";"
                                 + String(f_tl.i) + ";"
                                 + String(f_tl.it)+";"
                                 + String(f_tl.itt)+";"
                                 + String(f_tl.ittt)+";"
                                 + String(f_tl.r)+";"
                                 + String(f_tl.rd)+";"
                                 + String(f_tl.rdd)+";"
                                 + String(f_tl.rt)+";"
                                 + String(f_tl.rtt)+";"
                                 + String(f_tl.rtd)+";"
                                 + String(f_tl.rddd)+";"
                                 + String(f_tl.rtdd)+";"
                                 + String(f_tl.rttd)+";"
                                 + String(f_tl.rttt)+";",
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_tv.T) + ";"
                                 + String(f_tv.d) + ";"
                                 + String(f_tv.tau) + ";"
                                 + String(f_tv.delta) + ";"
                                 + String(f_tv.i) + ";"
                                 + String(f_tv.it)+";"
                                 + String(f_tv.itt)+";"
                                 + String(f_tv.ittt)+";"
                                 + String(f_tv.r)+";"
                                 + String(f_tv.rd)+";"
                                 + String(f_tv.rdd)+";"
                                 + String(f_tv.rt)+";"
                                 + String(f_tv.rtt)+";"
                                 + String(f_tv.rtd)+";"
                                 + String(f_tv.rddd)+";"
                                 + String(f_tv.rtdd)+";"
                                 + String(f_tv.rttd)+";"
                                 + String(f_tv.rttt)+";",
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_IIR.T) + ";"
                                 + String(f_IIR.d) + ";"
                                 + String(f_IIR.tau) + ";"
                                 + String(f_IIR.delta) + ";"
                                 + String(f_IIR.i) + ";"
                                 + String(f_IIR.it)+";"
                                 + String(f_IIR.itt)+";"
                                 + String(f_IIR.ittt)+";"
                                 + String(f_IIR.r)+";"
                                 + String(f_IIR.rd)+";"
                                 + String(f_IIR.rdd)+";"
                                 + String(f_IIR.rt)+";"
                                 + String(f_IIR.rtt)+";"
                                 + String(f_IIR.rtd)+";"
                                 + String(f_IIR.rddd)+";"
                                 + String(f_IIR.rtdd)+";"
                                 + String(f_IIR.rttd)+";"
                                 + String(f_IIR.rttt)+";",
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_ASHRAE.T) + ";"
                                 + String(f_ASHRAE.d) + ";"
                                 + String(f_ASHRAE.tau) + ";"
                                 + String(f_ASHRAE.delta) + ";"
                                 + String(f_ASHRAE.i) + ";"
                                 + String(f_ASHRAE.it)+";"
                                 + String(f_ASHRAE.itt)+";"
                                 + String(f_ASHRAE.ittt)+";"
                                 + String(f_ASHRAE.r)+";"
                                 + String(f_ASHRAE.rd)+";"
                                 + String(f_ASHRAE.rdd)+";"
                                 + String(f_ASHRAE.rt)+";"
                                 + String(f_ASHRAE.rtt)+";"
                                 + String(f_ASHRAE.rtd)+";"
                                 + String(f_ASHRAE.rddd)+";"
                                 + String(f_ASHRAE.rtdd)+";"
                                 + String(f_ASHRAE.rttd)+";"
                                 + String(f_ASHRAE.rttt)+";",
                                   fileName);
  Modelica.Utilities.Streams.print(String(f_NBP.T) + ";"
                                 + String(f_NBP.d) + ";"
                                 + String(f_NBP.tau) + ";"
                                 + String(f_NBP.delta) + ";"
                                 + String(f_NBP.i) + ";"
                                 + String(f_NBP.it)+";"
                                 + String(f_NBP.itt)+";"
                                 + String(f_NBP.ittt)+";"
                                 + String(f_NBP.r)+";"
                                 + String(f_NBP.rd)+";"
                                 + String(f_NBP.rdd)+";"
                                 + String(f_NBP.rt)+";"
                                 + String(f_NBP.rtt)+";"
                                 + String(f_NBP.rtd)+";"
                                 + String(f_NBP.rddd)+";"
                                 + String(f_NBP.rtdd)+";"
                                 + String(f_NBP.rttd)+";"
                                 + String(f_NBP.rttt)+";",
                                   fileName);
  Modelica.Utilities.Streams.print(String(f.T) + ";"
                                 + String(f.d) + ";"
                                 + String(f.tau) + ";"
                                 + String(f.delta) + ";"
                                 + String(f.i) + ";"
                                 + String(f.it)+";"
                                 + String(f.itt)+";"
                                 + String(f.ittt)+";"
                                 + String(f.r)+";"
                                 + String(f.rd)+";"
                                 + String(f.rdd)+";"
                                 + String(f.rt)+";"
                                 + String(f.rtt)+";"
                                 + String(f.rtd)+";"
                                 + String(f.rddd)+";"
                                 + String(f.rtdd)+";"
                                 + String(f.rttd)+";"
                                 + String(f.rttt)+";",
                                   fileName);

annotation (experiment(NumberOfIntervals=1));
end printHelmholtzDerivatives;
