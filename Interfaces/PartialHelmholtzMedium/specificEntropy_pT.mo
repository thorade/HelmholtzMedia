within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEntropy_pT
  "iteratively finds the specific entropy for a given p and T"

  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEntropy s "Specific Enthalpy";

protected
  MolarMass MM=fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta "reduced density";
  Real tau=T_crit/T "inverse reduced temperature";
  HelmholtzDerivs f;

  Density d;

algorithm
  assert(phase <> 2, "specificEntropy_pT error: pressure and temperature are not independent variables in two-phase state");
  d := density_pT(p=p, T=T, phase=phase);
  delta := d/d_crit;
  f.i := f_i(tau=tau, delta=delta);
  f.it := f_it(tau=tau, delta=delta);
  f.r := f_r(tau=tau, delta=delta);
  f.rt := f_rt(tau=tau, delta=delta);
  s := R*(tau*(f.it + f.rt) - f.i - f.r);

  // this is an iterative backward function
  // the two inverse functions are Temperature_ph and pressure_Th
  annotation (inverse(T=temperature_ps(p=p, s=s, phase=phase)));
end specificEntropy_pT;
