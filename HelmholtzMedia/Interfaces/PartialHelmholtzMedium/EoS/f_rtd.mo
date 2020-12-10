within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_rtd "residual part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_residual_delta_tau
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

  Real[nNonAna] Psi_t = {-a[i,10]*(2*tau-2)*Psi[i] for i in 1:nNonAna};
  Real[nNonAna] Theta_t = {-1 for i in 1:nNonAna};
  Real[nNonAna] Dis_t = {2*Theta[i]*Theta_t[i] for i in 1:nNonAna};
  Real[nNonAna] Disb_t = {if Dis[i]<>0 then a[i,6]*Dis[i]^(a[i,6]-1)*Dis_t[i] else 0 for i in 1:nNonAna};

  Real[nNonAna] Psi_d = {-a[i,9]*(2*delta-2)*Psi[i] for i in 1:nNonAna};
  Real[nNonAna] Theta_d = {if (delta-1)<>0 then (a[i,8]*(2*delta-2)*((delta-1)^2)^(0.5/a[i,7])) / (2*a[i,7]*(delta-1)^2) else 0 for i in 1:nNonAna};
  Real[nNonAna] Dis_d = {if (delta-1)<>0 then (2*a[i,11]*a[i,12]*((delta-1)^2)^a[i,12]) / (delta-1) +2*Theta[i]*Theta_d[i] else 0 for i in 1:nNonAna};
  Real[nNonAna] Disb_d = {if Dis[i]<>0 then a[i,6]*Dis[i]^(a[i,6]-1)*Dis_d[i] else 0 for i in 1:nNonAna};

  Real[nNonAna] Psi_td = {a[i,9]*a[i,10]*(2*delta-2)*(2*tau-2)*Psi[i] for i in 1:nNonAna};
  //Real[nNonAna] Theta_td = {0 for i in 1:nNonAna};
  Real[nNonAna] Dis_td = {2*Theta_d[i]*Theta_t[i] for i in 1:nNonAna};
  Real[nNonAna] Disb_td = {if Dis[i]<>0 then a[i,6]/Dis[i]^3* ((a[i,6]-1)*Dis[i]^(a[i,6]+1)*Dis_d[i]*Dis_t[i] + Dis[i]^(a[i,6]+2)*Dis_td[i]) else 0 for i in 1:nNonAna};

algorithm
  f_residual_delta_tau :=
      sum(p[i,1]*p[i,3]*p[i,2]*delta^(p[i,3] - 1)*tau^(p[i,2] - 1) for i in 1:nPoly)
    + sum(b[i,1]*b[i,2]*delta^(b[i,3] - 1)*tau^(b[i,2] - 1)*(b[i,3] - b[i,4]*delta^b[i,4])*exp(-delta^b[i,4]) for i in 1:nBwr)
    + sum(g[i,1]*delta^g[i,3]*tau^g[i,2]*exp(g[i,6]*(delta - g[i,9])^2 + g[i,7]*(tau - g[i,8])^2)*(g[i,3]/delta + 2*g[i,6]*(delta - g[i,9]))*(g[i,2]/tau + 2*g[i,7]*(tau - g[i,8])) for i in 1:nGauss)
    + sum(a[i,1]*(delta*Disb[i]*Psi_td[i] + delta*Psi[i]*Disb_td[i] + delta*Disb_d[i]*Psi_t[i] + delta*Disb_t[i]*Psi_d[i] + Disb[i]*Psi_t[i] + Psi[i]*Disb_t[i]) for i in 1:nNonAna);

end f_rtd;
