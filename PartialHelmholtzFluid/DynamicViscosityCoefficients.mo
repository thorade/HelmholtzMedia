within HelmholtzFluids.PartialHelmholtzFluid;
record DynamicViscosityCoefficients
  Temperature criticalTemperature "critical temperature";
  MolarVolume criticalMolarVolume "critical molar Volume";
  MolarMass molarMass "molar mass";

  // description of coeffs copied from RefProp
  constant Temperature epsilon_kappa "Lennard-Jones energy parameter";
  constant Real sigma "Lennard-Jones size parameter";
  constant Real[1,2] CET;

  constant Real[:,2] a "coefficients for S_mathfrak";
  constant Real[:,2] b "coefficients for second viscosity virial coefficent B";

  // coefficients for high density viscosity contribution
  constant Real[:,2] g "close-packed density delta_0";
  constant Real[:,5] e "simple poly";
  constant Real[:,5] nu_po "  numerator of rational poly";
  constant Real[:,5] de_po "denominator of rational poly";
  // constant Real[:,5] nu_ex "  numerator of exponential";
  // constant Real[:,5] de_ex "denominator of exponential";

end DynamicViscosityCoefficients;
