within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
record ThermalConductivityCoefficients
  // dilute gas / zero density terms
  constant Temperature reducingTemperature_0=1 "reducing temperature";
  constant Real reducingThermalConductivity_0=1 "usually unity";
  constant Real[:,2] lambda_0_coeffs "coeffs for dilute contribution";

  // residual / background terms
  constant Temperature reducingTemperature_residual=1 "reducing temperature";
  constant MolarVolume reducingMolarVolume_residual "reducing molar volume";
  constant Real reducingThermalConductivity_residual=1 "usually unity";
  constant Real[:,4] lambda_r_coeffs "coeffs for residual contribution";

  // critical enhancement terms
  constant Real nu=0.63 "universal exponent";
  constant Real gamma=1.239 "universal exponent";
  constant Real R0=1.03 "universal amplitude";
  constant Real z=0.063 "universal exponent, used for viscosity";
  constant Real c=1 "constant in viscosity, often set to 1";

  constant Real xi_0 "amplitude";
  constant Real Gamma_0 "amplitude";
  constant Real qd_inverse "modified effective cutoff parameter";
  constant Temperature T_ref "reference temperature";
end ThermalConductivityCoefficients;
