within HelmholtzMedia.Interfaces.Types;
record ThermalConductivityCoefficients
  constant HelmholtzMedia.Interfaces.Choices.ThermalConductivityModel thermalConductivityModel;
  constant HelmholtzMedia.Interfaces.Choices.ThermalConductivityCriticalEnhancementModel thermalConductivityCriticalEnhancementModel;
  // dilute gas / zero density terms
  constant Temperature reducingTemperature_0=1 "reducing temperature";
  constant Real reducingThermalConductivity_0=1 "usually unity";
  constant Real[:,2] lambda_0_num_coeffs = fill(0.0, 0, 2)
    "coeffs for dilute contribution numerator";
  constant Real[:,2] lambda_0_den_coeffs = fill(0.0, 0, 2)
    "coeffs for dilute contribution denominator";
  // background/residual terms
  constant Temperature reducingTemperature_background(min=1)=1
    "reducing temperature";
  constant MolarVolume reducingMolarVolume_background "reducing molar volume";
  constant Real reducingThermalConductivity_background=1 "usually unity";
  constant Real[:,4] lambda_b_coeffs = fill(0.0, 0, 4)
    "coeffs for background contribution";
  // critical enhancement terms
  constant Real nu=0.63 "universal exponent";
  constant Real gamma=1.239 "universal exponent";
  constant Real R0=1.03 "universal amplitude";
  constant Real z=0.063 "universal exponent, used for viscosity";
  constant Real c=1 "constant in viscosity, often set to 1";
  constant Real xi_0 "amplitude";
  constant Real Gamma_0 "amplitude";
  constant Real qd_inverse "modified effective cutoff parameter";
  constant Temperature T_ref(min=1,max=3e3) "reference temperature";
  // residual entropy scaling (RES)
  constant Real[:] RES_dilute = fill(0.0, 0)
    "RES dilute gas thermal conductivity in W/(m K): lambda_0 = RES_dilute[1]*T^4 + RES_dilute[2]*T^3 + ... + RES_dilute[5] (T in K)";
  constant Real[4] RES_n = zeros(4)
    "RES coefficients: lambda_r_plus = n1*x + n2*x^1.5 + n3*x^2 + n4*x^2.5 with x = s_plus/RES_xi";
  constant Real RES_xi = 1 "RES scaling factor of the reduced residual entropy";
end ThermalConductivityCoefficients;
