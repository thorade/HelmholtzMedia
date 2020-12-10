within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_rttt "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_residual_tau_tau
    "residual part of dimensionless Helmholtz energy";

protected
  constant Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  constant Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  constant Integer nGauss = size(helmholtzCoefficients.residualGauss,1);
  constant Integer nNonAna = size(helmholtzCoefficients.residualNonAnalytical,1);

  constant Real[nPoly,3] p = helmholtzCoefficients.residualPoly;
  constant Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  constant Real[nGauss,9] g = helmholtzCoefficients.residualGauss;
  constant Real[nNonAna,12] a = helmholtzCoefficients.residualNonAnalytical;

algorithm
  f_residual_tau_tau :=
      sum(delta^p[i,3]*p[i,1]*p[i,2]*tau^(p[i,2] - 3)*(p[i,2] - 1)*(p[i,2] - 2) for i in 1:nPoly)
    + sum(b[i,1]*b[i,2]*delta^b[i,3]*tau^(b[i,2] - 3)*exp(-delta^b[i,4])*(b[i,2] - 1)*(b[i,2] - 2) for i in 1:nBwr)
    + sum(delta^g[i,3]*g[i,1]*tau^(g[i,2] - 3)*exp(g[i,6]*(delta - g[i,9])^2 + g[i,7]*(g[i,8] - tau)^2)*(g[i,2]^3 - 6*g[i,2]^2*g[i,7]*g[i,8]*tau + 6*g[i,2]^2*g[i,7]*tau^2 - 3*g[i,2]^2 + 12*g[i,2]*g[i,7]^2*g[i,8]^2*tau^2 - 24*g[i,2]*g[i,7]^2*g[i,8]*tau^3 + 12*g[i,2]*g[i,7]^2*tau^4 + 6*g[i,2]*g[i,7]*g[i,8]*tau + 2*g[i,2] - 8*g[i,7]^3*g[i,8]^3*tau^3 + 24*g[i,7]^3*g[i,8]^2*tau^4 - 24*g[i,7]^3*g[i,8]*tau^5 + 8*g[i,7]^3*tau^6 - 12*g[i,7]^2*g[i,8]*tau^3 + 12*g[i,7]^2*tau^4) for i in 1:nGauss);

end f_rttt;
