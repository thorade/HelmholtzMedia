within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_r "residual part of dimensionless Helmholtz energy"

  input Real delta(min=0);
  input Real tau(min=0);
  output Real f_residual "residual part of dimensionless Helmholtz energy";

protected
  constant Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  constant Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  constant Integer nGauss = size(helmholtzCoefficients.residualGauss,1);
  constant Integer nNonAna = size(helmholtzCoefficients.residualNonAnalytical,1);

  constant Real[nPoly,3] p = helmholtzCoefficients.residualPoly;
  constant Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  constant Real[nGauss,9] g = helmholtzCoefficients.residualGauss;
  constant Real[nNonAna,12] a = helmholtzCoefficients.residualNonAnalytical;

  Real[nNonAna] Psi = {exp(-a[i,9]*(delta-1)^2 -a[i,10]*(tau-1)^2) for i in 1:nNonAna} "Psi";
  Real[nNonAna] Theta = {(1-tau) + a[i,8]*((delta-1)^2)^(0.5/a[i,7]) for i in 1:nNonAna} "Theta";
  Real[nNonAna] Dis = {Theta[i]^2 +a[i,11]*((delta-1)^2)^a[i,12] for i in 1:nNonAna} "Distance function";
  Real[nNonAna] Disb = {Dis[i]^a[i,6] for i in 1:nNonAna} "Distance function to the power of b";

algorithm
  f_residual :=
    sum(p[i,1]*tau^p[i,2]*delta^p[i,3] for i in 1:nPoly)
  + sum(b[i,1]*tau^b[i,2]*delta^b[i,3]*exp(-delta^b[i,4]) for i in 1:nBwr)
  + sum(g[i,1]*tau^g[i,2]*delta^g[i,3]*exp(g[i,6]*(delta - g[i,9])^g[i,4] + g[i,7]*(tau - g[i,8])^g[i,5]) for i in 1:nGauss)
  + sum(a[i,1]*Disb[i]*delta*Psi[i] for i in 1:nNonAna);
end f_r;
