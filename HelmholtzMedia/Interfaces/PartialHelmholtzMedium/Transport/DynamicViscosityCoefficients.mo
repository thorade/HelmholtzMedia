within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
record DynamicViscosityCoefficients
  constant DynamicViscosityModel dynamicViscosityModel;
  constant CollisionIntegralModel collisionIntegralModel;
  // collision integral S_mathfrak
  constant Temperature epsilon_kappa(min=0) "Lennard-Jones energy parameter";
  constant Real sigma "Lennard-Jones size parameter";
  constant Real[:,2] a = fill(0.0, 0, 2) "coefficients for collision integral";
  // zero density dependence, via Chapman-Enskog-Theory
  constant Real[:,2] CET = fill(0.0, 0, 2);
  constant Temperature reducingTemperature_0(min=1)=1 "reducing temperature";
  constant Real reducingViscosity_0=1 "reducing viscosity";
  // initial density dependence, via Rainwater-Friend-Theory
  constant Temperature reducingTemperature_1(min=1)=1 "reducing temperature";
  constant Real reducingViscosity_1=1 "reducing viscosity";
  constant Real[:,2] b = fill(0.0, 0, 2)
    "coefficients for second viscosity virial coefficent B";
  // residual / high density viscosity contribution
  constant Real[:,1] c =  fill(0.0, 0, 1)
    "coefficients for residual viscosity contribution in VS2 model";
  constant Temperature reducingTemperature_residual(min=1)=1
    "reducing temperature";
  constant MolarVolume reducingMolarVolume_residual=1 "reducing molar volume";
  constant Real reducingViscosity_residual=1 "reducing viscosity";
  constant Real[:,2] g = fill(0.0, 0, 2) "close-packed density delta_0";
  constant Real[:,5] e = fill(0.0, 0, 5) "simple poly";
  constant Real[:,5] nu_po = fill(0.0, 0, 5) "  numerator of rational poly";
  constant Real[:,5] de_po = fill(0.0, 0, 5) "denominator of rational poly";
  // constant Real[:,5] nu_ex = fill(0.0, 0, 5) "  numerator of exponential";
  // constant Real[:,5] de_ex = fill(0.0, 0, 5) "denominator of exponential";
end DynamicViscosityCoefficients;
