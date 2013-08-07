within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_r "residual part of dimensionless Helmholtz energy"

  input Real delta(min=0);
  input Real tau(min=0);
  output Real f_residual "residual part of dimensionless Helmholtz energy";

protected
  final constant Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  final constant Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  final constant Integer nGauss = size(helmholtzCoefficients.residualGauss,1);

  final constant Real[nPoly,4] p = helmholtzCoefficients.residualPoly;
  final constant Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  final constant Real[nGauss,9] g = helmholtzCoefficients.residualGauss;

algorithm
  f_residual :=
    sum(p[i,1]*tau^p[i,2]*delta^p[i,3] for i in 1:nPoly)
  + sum(b[i,1]*tau^b[i,2]*delta^b[i,3]*exp(-delta^b[i,4]) for i in 1:nBwr)
  + sum(g[i,1]*tau^g[i,2]*delta^g[i,3]*exp(g[i,6]*(delta - g[i,9])^g[i,4] + g[i,7]*(tau - g[i,8])^g[i,5]) for i in 1:nGauss);
end f_r;
