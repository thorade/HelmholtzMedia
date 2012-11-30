within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_rddd "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_residual_delta_delta_delta
    "residual part of dimensionless Helmholtz energy";

protected
  Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  Integer nGauss = size(helmholtzCoefficients.residualGauss,1);

  Real[nPoly,4] p = helmholtzCoefficients.residualPoly;
  Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  Real[nGauss,9] g = helmholtzCoefficients.residualGauss;

algorithm
  f_residual_delta_delta_delta :=
      sum(delta^(p[i,3] - 3)*p[i,1]*p[i,3]*tau^p[i,2]*(p[i,3] - 1)*(p[i,3] - 2) for i in 1:nPoly)
    + sum(-b[i,1]*delta^(b[i,3] - 3)*tau^b[i,2]*exp(-delta^b[i,4])*(b[i,4]^3*delta^b[i,4] - 3*b[i,4]^2*delta^b[i,4] - 2*b[i,3] + 3*b[i,4]^2*delta^(2*b[i,4]) - 3*b[i,4]^3*delta^(2*b[i,4]) + b[i,4]^3*delta^(3*b[i,4]) + 3*b[i,3]^2 - b[i,3]^3 + 2*b[i,4]*delta^b[i,4] - 3*b[i,3]*b[i,4]^2*delta^(2*b[i,4]) - 6*b[i,3]*b[i,4]*delta^b[i,4] + 3*b[i,3]*b[i,4]^2*delta^b[i,4] + 3*b[i,3]^2*b[i,4]*delta^b[i,4]) for i in 1:nBwr)
    + sum(delta^(g[i,3] - 3)*g[i,1]*tau^g[i,2]*exp(g[i,6]*(delta - g[i,9])^2 + g[i,7]*(g[i,8] - tau)^2)*(8*delta^6*g[i,6]^3 - 24*delta^5*g[i,6]^3*g[i,9] + 12*delta^4*g[i,3]*g[i,6]^2 + 24*delta^4*g[i,6]^3*g[i,9]^2 + 12*delta^4*g[i,6]^2 - 24*delta^3*g[i,3]*g[i,6]^2*g[i,9] - 8*delta^3*g[i,6]^3*g[i,9]^3 - 12*delta^3*g[i,6]^2*g[i,9] + 6*delta^2*g[i,3]^2*g[i,6] + 12*delta^2*g[i,3]*g[i,6]^2*g[i,9]^2 - 6*delta*g[i,3]^2*g[i,6]*g[i,9] + 6*delta*g[i,3]*g[i,6]*g[i,9] + g[i,3]^3 - 3*g[i,3]^2 + 2*g[i,3]) for i in 1:nGauss);

end f_rddd;
