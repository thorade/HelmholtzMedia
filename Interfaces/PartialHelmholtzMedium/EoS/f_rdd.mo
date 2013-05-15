within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_rdd "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_residual_delta_delta
    "residual part of dimensionless Helmholtz energy";

protected
  final constant Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  final constant Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  final constant Integer nGauss = size(helmholtzCoefficients.residualGauss,1);

  final constant Real[nPoly,4] p = helmholtzCoefficients.residualPoly;
  final constant Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  final constant Real[nGauss,9] g = helmholtzCoefficients.residualGauss;

algorithm
  f_residual_delta_delta :=
      sum(p[i,1]*p[i,3]*(p[i,3] - 1)*delta^(p[i,3] - 2)*tau^p[i,2] for i in 1:nPoly)
    + sum(b[i,1]*exp(-delta^b[i,4])*(delta^(b[i,3] - 2)*tau^b[i,2]*((b[i,3] - b[i,4]*delta^b[i,4])*(b[i,3] - 1 - b[i,4]*delta^b[i,4]) - b[i,4]^2*delta^b[i,4])) for i in 1:nBwr)
    + sum(g[i,1]*tau^g[i,2]*exp(g[i,6]*(delta - g[i,9])^2 + g[i,7]*(tau - g[i,8])^2)*(  2*g[i,6]*delta^g[i,3] + 4*g[i,6]^2*delta^g[i,3]*(delta - g[i,9])^2 + 4*g[i,3]*g[i,6]*delta^(g[i,3] - 1)*(delta - g[i,9]) +   g[i,3]*(g[i,3] - 1)*delta^(g[i,3] - 2)) for i in 1:nGauss);

end f_rdd;
