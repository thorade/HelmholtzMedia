within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
record DynamicViscosityCoefficients

  constant Types.DynamicViscosityModel dynamicViscosityModel;

  // description of coeffs copied from RefProp
  constant Temperature epsilon_kappa "Lennard-Jones energy parameter";
  constant Real sigma "Lennard-Jones size parameter";
  constant Real[1,2] CET;

  // zero density dependence /collision integral S_mathfrak
  constant Real[:,2] a "coefficients for collision integral";

  // initial density dependence
  constant Temperature reducingTemperature_0=1 "reducing temperature";
  constant Real reducingViscosity_0=1 "reducing viscosity";
  constant Real[:,2] b "coefficients for second viscosity virial coefficent B";

  // residual / high density viscosity contribution
  constant Temperature reducingTemperature_residual "reducing temperature";
  constant MolarVolume reducingMolarVolume_residual "reducing molar volume";
  constant Real reducingViscosity_residual "reducing viscosity";
  constant Boolean hasGeneralizedDelta0=true;
  constant Real[:,2] g "close-packed density delta_0";
  constant Real[:,5] e "simple poly";
  constant Real[:,5] nu_po "  numerator of rational poly";
  constant Real[:,5] de_po "denominator of rational poly";
  // constant Real[:,5] nu_ex "  numerator of exponential";
  // constant Real[:,5] de_ex "denominator of exponential";
end DynamicViscosityCoefficients;
