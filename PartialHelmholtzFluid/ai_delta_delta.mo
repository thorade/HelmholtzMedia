within HelmholtzFluids.PartialHelmholtzFluid;
function ai_delta_delta "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real alpha_ideal_delta_delta
    "ideal part of dimensionless Helmholtz energy";

algorithm
    alpha_ideal_delta_delta := -1.0/(delta*delta);
end ai_delta_delta;
