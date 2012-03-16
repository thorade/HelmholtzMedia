within HelmholtzMedia.Examples;
model Validate_Derivatives_SinglePhase
  "compare analytical derivatives to numerical derivatives"

  package medium = HelmholtzFluids.R134a;
  // choose d and T which will result in single-phase
  parameter medium.Density d=2e+3;
  parameter medium.Temperature T=298.15;

  medium.ThermodynamicState state;
  medium.HelmholtzDerivs f;
  medium.ThermodynamicState d_plus;
  medium.ThermodynamicState d_minus;
  medium.ThermodynamicState T_plus;
  medium.ThermodynamicState T_minus;

  // Enthalpy derivatives
  medium.Types.DerEnthalpyByDensity dhdT_analytical;
  medium.Types.DerEnthalpyByDensity dhdT_numerical;
  medium.Types.DerEnthalpyByTemperature dhTd_analytical;
  medium.Types.DerEnthalpyByTemperature dhTd_numerical;
  // Energy derivatives
  medium.Types.DerEnergyByDensity dudT_analytical;
  medium.Types.DerEnergyByDensity dudT_numerical;
  medium.Types.DerEnergyByTemperature duTd_analytical;
  medium.Types.DerEnergyByTemperature duTd_numerical;
  // Entropy derivatives
  medium.Types.DerEntropyByDensity dsdT_analytical;
  medium.Types.DerEntropyByDensity dsdT_numerical;
  medium.Types.DerEntropyByTemperature dsTd_analytical;
  medium.Types.DerEntropyByTemperature dsTd_numerical;
  // Gibbs derivatives
  medium.Types.DerEnergyByDensity dgdT_analytical;
  medium.Types.DerEnergyByDensity dgdT_numerical;
  medium.Types.DerEnergyByTemperature dgTd_analytical;
  medium.Types.DerEnergyByTemperature dgTd_numerical;

equation
  state=medium.setState_dTX(d=d, T=T);
  f=medium.setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
  d_plus=medium.setState_dTX(d=d*1.00001, T=T);
  d_minus=medium.setState_dTX(d=d*0.99999, T=T);
  T_plus=medium.setState_dTX(d=d, T=T*1.00001);
  T_minus=medium.setState_dTX(d=d, T=T*0.99999);

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|");

  Modelica.Utilities.Streams.print("Enthalpy");
  // check (dh/dd)@T=const
  dhdT_analytical = medium.specificEnthalpy_derd_T(state=state, f=f);
  dhdT_numerical = (d_plus.h-d_minus.h)/(d_plus.d-d_minus.d);
  Modelica.Utilities.Streams.print("(dh/dd)@T=const analytical= " + String(dhdT_analytical));
  Modelica.Utilities.Streams.print("(dh/dd)@T=const  numerical= " + String(dhdT_numerical));
  // check (dh/dT)@d=const
  dhTd_analytical = medium.specificEnthalpy_derT_d(state=state, f=f);
  dhTd_numerical = (T_plus.h-T_minus.h)/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("(dh/dT)@d=const analytical= " + String(dhTd_analytical));
  Modelica.Utilities.Streams.print("(dh/dT)@d=const  numerical= " + String(dhTd_numerical));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("internal Energy");
  // check (du/dd)@T=const
  dudT_analytical = f.R*T/d*f.tau*f.delta*f.rtd;
  dudT_numerical = (d_plus.u-d_minus.u)/(d_plus.d-d_minus.d);
  Modelica.Utilities.Streams.print("(du/dd)@T=const analytical= " + String(dudT_analytical));
  Modelica.Utilities.Streams.print("(du/dd)@T=const  numerical= " + String(dudT_numerical));
  // check (du/dT)@d=const
  duTd_analytical = medium.specificHeatCapacityCv(state=state, f=f);
  duTd_numerical = (T_plus.u-T_minus.u)/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("(du/dT)@d=const analytical= " + String(duTd_analytical));
  Modelica.Utilities.Streams.print("(du/dT)@d=const  numerical= " + String(duTd_numerical));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Entropy");
  // check (ds/dd)@T=const
  dsdT_analytical = f.R/d*(-(1+f.delta*f.rd)+(0+f.tau*f.delta*f.rtd));
  dsdT_numerical = (d_plus.s-d_minus.s)/(d_plus.d-d_minus.d);
  Modelica.Utilities.Streams.print("(ds/dd)@T=const analytical= " + String(dsdT_analytical));
  Modelica.Utilities.Streams.print("(ds/dd)@T=const  numerical= " + String(dsdT_numerical));
  // check (ds/dT)@d=const
  dsTd_analytical = f.R/T*(-f.tau^2*(f.itt+f.rtt));
  dsTd_numerical = (T_plus.s-T_minus.s)/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("(ds/dT)@d=const analytical= " + String(dsTd_analytical));
  Modelica.Utilities.Streams.print("(ds/dT)@d=const  numerical= " + String(dsTd_numerical));

  Modelica.Utilities.Streams.print(" ");
  Modelica.Utilities.Streams.print("Gibbs energy");
  // check (dg/dd)@T=const
  dgdT_analytical = f.R*T/d*(1+2*f.delta*f.rd + f.delta^2*f.rdd);
  dgdT_numerical = ((d_plus.h-d_plus.T*d_plus.s)-(d_minus.h-d_minus.T*d_minus.s))/(d_plus.d-d_minus.d);
  Modelica.Utilities.Streams.print("(dg/dd)@T=const analytical= " + String(dgdT_analytical));
  Modelica.Utilities.Streams.print("(dg/dd)@T=const  numerical= " + String(dgdT_numerical));
  // check (dg/dT)@d=const
  dgTd_analytical = f.R*(f.i+f.r+1+f.delta*f.rd -f.tau*(f.it+f.rt) - f.tau*f.delta*f.rtd);
  dgTd_numerical = ((T_plus.h-T_plus.T*T_plus.s)-(T_minus.h-T_minus.T*T_minus.s))/(T_plus.T-T_minus.T);
  Modelica.Utilities.Streams.print("(dg/dT)@d=const analytical= " + String(dgTd_analytical));
  Modelica.Utilities.Streams.print("(dg/dT)@d=const  numerical= " + String(dgTd_numerical));

  annotation (experiment(NumberOfIntervals=1), __Dymola_experimentSetupOutput);
end Validate_Derivatives_SinglePhase;
