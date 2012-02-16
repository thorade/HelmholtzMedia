within HelmholtzFluids.PartialHelmholtzFluid;
record DynamicViscosityCoefficients

  // description of coeffs copied from RefProp
  constant Temperature epsilon_kappa "Lennard-Jones energy parameter";
  constant Real sigma "Lennard-Jones size parameter";
  constant Real[1,2] CET;

  // zero density dependence /collision integral S_mathfrak
  constant Real[:,2] a "coefficients for collision integral";

  // initial density dependence
  constant Real[:,2] b "coefficients for second viscosity virial coefficent B";

  // high density viscosity contribution
  Temperature reducingTemperature "reducing temperature";
  MolarVolume reducingMolarVolume "reducing molar volume";
  constant Real[:,2] g "close-packed density delta_0";
  constant Real[:,5] e "simple poly";
  constant Real[:,5] nu_po "  numerator of rational poly";
  constant Real[:,5] de_po "denominator of rational poly";
  // constant Real[:,5] nu_ex "  numerator of exponential";
  // constant Real[:,5] de_ex "denominator of exponential";

end DynamicViscosityCoefficients;
