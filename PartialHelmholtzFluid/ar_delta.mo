within HelmholtzFluids.PartialHelmholtzFluid;
function ar_delta "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real alpha_residual_delta
    "residual part of dimensionless Helmholtz energy";

protected
  Real[size(helmholtzCoefficients.n_residual,1)] n=helmholtzCoefficients.n_residual;
  Real[size(helmholtzCoefficients.c,1)] c=helmholtzCoefficients.c;
  Real[size(helmholtzCoefficients.d,1)] d=helmholtzCoefficients.d;
  Real[size(helmholtzCoefficients.t,1)] t=helmholtzCoefficients.t;
  Real[size(helmholtzCoefficients.crit_epsilon,1)] epsilon=helmholtzCoefficients.crit_epsilon;
  Real[size(helmholtzCoefficients.crit_beta,1)] beta=helmholtzCoefficients.crit_beta;
  Real[size(helmholtzCoefficients.crit_eta,1)] eta=helmholtzCoefficients.crit_eta;
  Real[size(helmholtzCoefficients.crit_gamma,1)] gamma=helmholtzCoefficients.crit_gamma;

algorithm
  alpha_residual_delta := sum(n[i]*d[i]*delta^(d[i] - 1)*tau^t[i] for i in 1:7)
    + sum(n[i]*exp(-delta^c[i])*(delta^(d[i] - 1)*tau^t[i]*(d[i] - c[i]*delta^c[i])) for i in 8:23)
    + sum(n[i]*delta^d[i]*tau^t[i]*exp(-eta[i]*(delta - epsilon[i])^2 - beta[i]*(tau - gamma[i])^2)*(d[i]/delta - 2*eta[i]*(delta - epsilon[i])) for i in 24:25);

end ar_delta;
