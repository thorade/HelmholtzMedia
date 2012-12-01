within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_rttd "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_residual_tau_tau_delta
    "residual part of dimensionless Helmholtz energy";

protected
  Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  Integer nGauss = size(helmholtzCoefficients.residualGauss,1);

  Real[nPoly,4] p = helmholtzCoefficients.residualPoly;
  Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  Real[nGauss,9] g = helmholtzCoefficients.residualGauss;

algorithm
  f_residual_tau_tau_delta :=
      sum(delta^(p[i,3] - 1)*p[i,1]*p[i,2]*p[i,3]*tau^(p[i,2] - 2)*(p[i,2] - 1) for i in 1:nPoly)
    + sum(b[i,1]*b[i,2]*delta^(b[i,3] - 1)*tau^(b[i,2] - 2)*exp(-delta^b[i,4])*(b[i,3] - b[i,4]*delta^b[i,4])*(b[i,2] - 1) for i in 1:nBwr)
    + sum(delta^(g[i,3] - 1)*g[i,1]*tau^(g[i,2] - 2)*exp(g[i,6]*(delta - g[i,9])^2 + g[i,7]*(g[i,8] - tau)^2)*(2*g[i,6]*delta^2 - 2*g[i,6]*g[i,9]*delta + g[i,3])*(g[i,2]^2 - 4*g[i,2]*g[i,7]*g[i,8]*tau + 4*g[i,2]*g[i,7]*tau^2 - g[i,2] + 4*g[i,7]^2*g[i,8]^2*tau^2 - 8*g[i,7]^2*g[i,8]*tau^3 + 4*g[i,7]^2*tau^4 + 2*g[i,7]*tau^2) for i in 1:nGauss);

end f_rttd;
