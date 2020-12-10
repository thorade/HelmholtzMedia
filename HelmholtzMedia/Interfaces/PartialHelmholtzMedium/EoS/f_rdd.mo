within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_rdd "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_residual_delta_delta
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

  Real[nNonAna] Psi = {exp(-a[i,9]*(delta-1)^2 -a[i,10]*(tau-1)^2) for i in 1:nNonAna} "Psi";
  Real[nNonAna] Theta = {(1-tau) + a[i,8]*((delta-1)^2)^(0.5/a[i,7]) for i in 1:nNonAna} "Theta";
  Real[nNonAna] Dis = {Theta[i]^2 +a[i,11]*((delta-1)^2)^a[i,12] for i in 1:nNonAna} "Distance function";
  Real[nNonAna] Disb = {Dis[i]^a[i,6] for i in 1:nNonAna} "Distance function to the power of b";

  Real[nNonAna] Psi_d = {-a[i,9]*(2*delta-2)*Psi[i] for i in 1:nNonAna};
  Real[nNonAna] Theta_d = {if (delta-1)<>0 then (a[i,8]*(2*delta-2)*((delta-1)^2)^(0.5/a[i,7])) / (2*a[i,7]*(delta-1)^2) else 0 for i in 1:nNonAna};
  Real[nNonAna] Dis_d = {if (delta-1)<>0 then (2*a[i,11]*a[i,12]*((delta-1)^2)^a[i,12]) / (delta-1) +2*Theta[i]*Theta_d[i] else 0 for i in 1:nNonAna};
  Real[nNonAna] Disb_d = {if Dis[i]<>0 then a[i,6]*Dis[i]^(a[i,6]-1)*Dis_d[i] else 0 for i in 1:nNonAna};

  Real[nNonAna] Psi_dd = {2*a[i,9]*(2*a[i,9]*(delta-1)^2-1)*Psi[i] for i in 1:nNonAna};
  Real[nNonAna] Theta_dd = {if (delta-1)<>0 then (a[i,8]*(1/a[i,7]-1)*((delta-1)^2)^(0.5/a[i,7])) / (a[i,7]*(delta-1)^2) else 0 for i in 1:nNonAna};
  Real[nNonAna] Dis_dd = {if (delta-1)<>0 then 2*Theta[i]*Theta_dd[i] + 2*Theta_d[i]^2 + (2*((delta-1)^2)^a[i,12]) / (delta-1)^2 *(2*a[i,11]*a[i,12]^2-a[i,11]*a[i,12]) else 0 for i in 1:nNonAna};
  Real[nNonAna] Disb_dd = {if (delta-1)<>0 then a[i,6]*(a[i,6]*Dis_d[i]^2 + Dis[i]*Dis_dd[i] - Dis_d[i]^2)*Dis[i]^(a[i,6]-2) else 0 for i in 1:nNonAna};

algorithm
  f_residual_delta_delta :=
      sum(p[i,1]*p[i,3]*(p[i,3] - 1)*delta^(p[i,3] - 2)*tau^p[i,2] for i in 1:nPoly)
    + sum(b[i,1]*exp(-delta^b[i,4])*(delta^(b[i,3] - 2)*tau^b[i,2]*((b[i,3] - b[i,4]*delta^b[i,4])*(b[i,3] - 1 - b[i,4]*delta^b[i,4]) - b[i,4]^2*delta^b[i,4])) for i in 1:nBwr)
    + sum(g[i,1]*tau^g[i,2]*exp(g[i,6]*(delta - g[i,9])^2 + g[i,7]*(tau - g[i,8])^2)*(  2*g[i,6]*delta^g[i,3] + 4*g[i,6]^2*delta^g[i,3]*(delta - g[i,9])^2 + 4*g[i,3]*g[i,6]*delta^(g[i,3] - 1)*(delta - g[i,9]) +   g[i,3]*(g[i,3] - 1)*delta^(g[i,3] - 2)) for i in 1:nGauss)
    + sum(a[i,1]*(delta*Disb[i]*Psi_dd[i] + delta*Psi[i]*Disb_dd[i] + 2*delta*Disb_d[i]*Psi_d[i] + 2*Disb[i]*Psi_d[i] + 2*Psi[i]*Disb_d[i]) for i in 1:nNonAna);

end f_rdd;
