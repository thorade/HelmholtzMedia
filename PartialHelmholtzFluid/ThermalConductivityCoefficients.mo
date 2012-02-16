within HelmholtzFluids.PartialHelmholtzFluid;
record ThermalConductivityCoefficients
  Temperature reducingTemperature
    "reducing temperature (either 1 or some value close to critical temperature)";
  MolarVolume reducingMolarVolume "reducing molar volume";
  Real reducingThermalConductivity=1 "usually unity, sometimes different value";
  constant Real[:,2] lambda_0_coeffs "coeffs for dilute contribution";
  constant Real[:,4] lambda_r_coeffs "coeffs for residual contribution";

  Real nu "universal exponent";
  Real gamma "universal exponent";
  Real R0 "universal amplitude";
  Real z "universal exponent, used for viscosity";
  Real c "constant in viscosity, often set to 1";
  Real xi_0 "amplitude";
  Real Gamma_0 "amplitude";
  Real qd_inverse "modified effective cutoff parameter";
  Temperature T_ref "reference temperature";

end ThermalConductivityCoefficients;
