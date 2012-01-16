within HelmholtzFluids.PartialHelmholtzFluid;
record HelmholtzCoefficients
  "Coefficients for Helmholtz energy equations of state"

  //ideal gas: n-butane substance specific coefficients
  constant Real[:] n_ideal;
  constant Real[:] Theta;

  //residual part: n-butane substance specific coefficients
  constant Real[:] n_residual;

  //residual part functional form exponents
  constant Integer[:] c "dimensionless density exponents in exponential term";
  constant Integer[:] d "dimensionless density exponents";
  constant Real[:] t "dimensionless temperature exponents";

  constant Integer[:] crit_eta
    "dimensionless coefficients for critical region correction";
  constant Integer[:] crit_beta
    "dimensionless coefficients for critical region correction";
  constant Real[:] crit_epsilon
    "dimensionless coefficients for critical region correction";
  constant Real[:] crit_gamma
    "dimensionless coefficients for critical region correction";
end HelmholtzCoefficients;
