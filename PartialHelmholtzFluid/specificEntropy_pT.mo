within HelmholtzFluids.PartialHelmholtzFluid;
function specificEntropy_pT
  "iteratively finds the specific entropy for a given p and T"

  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEntropy s "Specific Enthalpy";

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta "reduced density";
  Real tau=T_crit/T "inverse reduced temperature";

  Density d;

algorithm
  assert(phase <> 2, "specificEntropy_pT error: pressure and temperature are not independent variables in two-phase state");
  d := density_pT(p=p, T=T, phase=phase);
  delta := d/d_crit;
  s := (tau*(ai_tau(delta=delta, tau=tau) + ar_tau(delta=delta, tau=tau)) - ai(delta=delta, tau=tau) - ar(delta=delta, tau=tau))*R;

  // this is an iterative backward function
  // pressure_dT is the corresponding forward function
  // annotation (inverse(p=pressure_dT(d=d, T=T, phase=phase)));
end specificEntropy_pT;
