within HelmholtzFluids.PartialHelmholtzFluid;
function ai_tau "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real alpha_ideal_tau "ideal part of dimensionless Helmholtz energy";

protected
  Real[size(helmholtzCoefficients.n_ideal,1)] n=helmholtzCoefficients.n_ideal;
  Real[size(helmholtzCoefficients.Theta,1)] Theta=helmholtzCoefficients.Theta;

algorithm
  alpha_ideal_tau := n[2] + n[3]/tau
    + sum(n[i]*Theta[i]*((1 - exp(-Theta[i]*tau))^(-1) - 1) for i in 4:7);
end ai_tau;
