within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_r_derivs
  "residual part of the dimensionless Helmholtz energy and its derivatives up to the given order, in one pass over all terms"
  input Real delta "reduced density";
  input Real tau "inverse reduced temperature";
  input Integer order(min=0, max=3)=2 "highest order of the derivatives to compute";
  output Real r "f_r";
  output Real rd "(df_r/ddelta)@tau=const";
  output Real rdd "(d2f_r/ddelta2)@tau=const";
  output Real rt "(df_r/dtau)@delta=const";
  output Real rtt "(d2f_r/dtau2)@delta=const";
  output Real rtd "(d2f_r/dtau ddelta)";
  output Real rttt "(d3f_r/dtau3)@delta=const";
  output Real rttd "(d3f_r/dtau2 ddelta)";
  output Real rtdd "(d3f_r/dtau ddelta2)";
  output Real rddd "(d3f_r/ddelta3)@tau=const";

protected
  constant Integer nPoly = size(helmholtzCoefficients.residualPoly,1);
  constant Integer nBwr = size(helmholtzCoefficients.residualBwr,1);
  constant Integer nGauss = size(helmholtzCoefficients.residualGauss,1);
  constant Integer nNonAna = size(helmholtzCoefficients.residualNonAnalytical,1);

  constant Real[nPoly,3] p = helmholtzCoefficients.residualPoly;
  constant Real[nBwr,4] b = helmholtzCoefficients.residualBwr;
  constant Real[nGauss,9] g = helmholtzCoefficients.residualGauss;
  constant Real[nNonAna,12] a = helmholtzCoefficients.residualNonAnalytical;

  Real del "reduced density, protected against zero";
  Real lnd "log(delta)";
  Real lnt "log(tau)";
  Real k "value of the current term";
  Real D1 "delta*d(term)/ddelta / term";
  Real D2 "delta^2*d2(term)/ddelta2 / term";
  Real D3 "delta^3*d3(term)/ddelta3 / term";
  Real T1 "tau*d(term)/dtau / term";
  Real T2 "tau^2*d2(term)/dtau2 / term";
  Real T3 "tau^3*d3(term)/dtau3 / term";
  Real x "c*delta^c (exponential terms)";
  Real Rd "2*eta*(2*delta^2-eps*delta) (Gaussian terms)";
  Real Rt "2*beta*(2*tau^2-gamma*tau) (Gaussian terms)";

  // non-analytical terms (only first and second derivatives, as in f_rd ... f_rtd)
  Real Psi;
  Real Theta;
  Real Dis;
  Real Disb;
  Real Psi_d;
  Real Theta_d;
  Real Dis_d;
  Real Disb_d;
  Real Psi_dd;
  Real Theta_dd;
  Real Dis_dd;
  Real Disb_dd;
  Real Psi_t;
  Real Dis_t;
  Real Disb_t;
  Real Psi_tt;
  Real Disb_tt;
  Real Psi_td;
  Real Dis_td;
  Real Disb_td;
  Real na_r=0;
  Real na_rd=0;
  Real na_rdd=0;
  Real na_rt=0;
  Real na_rtt=0;
  Real na_rtd=0;

algorithm
  // The sums are accumulated as delta^m*tau^n*(d^(m+n) f_r / ddelta^m dtau^n), which for each term is
  // the term value k times factors D_m (delta part) and T_n (tau part); the powers delta^d*tau^t are computed
  // once per term as exp(d*log(delta) + t*log(tau)). The division by delta^m*tau^n follows at the end.
  r := 0;
  rd := 0;
  rdd := 0;
  rt := 0;
  rtt := 0;
  rtd := 0;
  rttt := 0;
  rttd := 0;
  rtdd := 0;
  rddd := 0;
  del := max(delta, Modelica.Constants.small);
  lnd := log(del);
  lnt := log(tau);

  // polynomial terms n*delta^d*tau^t
  for i in 1:nPoly loop
    k := p[i,1]*exp(p[i,3]*lnd + p[i,2]*lnt);
    r := r + k;
    if order >= 1 then
      D1 := p[i,3];
      T1 := p[i,2];
      rd := rd + k*D1;
      rt := rt + k*T1;
      if order >= 2 then
        D2 := D1*(D1 - 1);
        T2 := T1*(T1 - 1);
        rdd := rdd + k*D2;
        rtt := rtt + k*T2;
        rtd := rtd + k*D1*T1;
        if order >= 3 then
          D3 := D2*(D1 - 2);
          T3 := T2*(T1 - 2);
          rddd := rddd + k*D3;
          rttt := rttt + k*T3;
          rttd := rttd + k*T2*D1;
          rtdd := rtdd + k*D2*T1;
        end if;
      end if;
    end if;
  end for;

  // exponential terms n*delta^d*tau^t*exp(-delta^c)
  for i in 1:nBwr loop
    x := exp(b[i,4]*lnd) "delta^c";
    k := b[i,1]*exp(b[i,3]*lnd + b[i,2]*lnt - x);
    x := b[i,4]*x "c*delta^c";
    r := r + k;
    if order >= 1 then
      D1 := b[i,3] - x;
      T1 := b[i,2];
      rd := rd + k*D1;
      rt := rt + k*T1;
      if order >= 2 then
        D2 := D1*(D1 - 1) - b[i,4]*x;
        T2 := T1*(T1 - 1);
        rdd := rdd + k*D2;
        rtt := rtt + k*T2;
        rtd := rtd + k*D1*T1;
        if order >= 3 then
          D3 := D1*D2 - 2*D1*b[i,4]*x + b[i,4]*x - b[i,4]^2*x - 2*D2;
          T3 := T2*(T1 - 2);
          rddd := rddd + k*D3;
          rttt := rttt + k*T3;
          rttd := rttd + k*T2*D1;
          rtdd := rtdd + k*D2*T1;
        end if;
      end if;
    end if;
  end for;

  // Gaussian bell-shaped terms n*delta^d*tau^t*exp(eta*(delta-eps)^2 + beta*(tau-gamma)^2)
  for i in 1:nGauss loop
    k := g[i,1]*exp(g[i,3]*lnd + g[i,2]*lnt + g[i,6]*(del - g[i,9])^2 + g[i,7]*(tau - g[i,8])^2);
    r := r + k;
    if order >= 1 then
      D1 := g[i,3] + 2*g[i,6]*del*(del - g[i,9]);
      T1 := g[i,2] + 2*g[i,7]*tau*(tau - g[i,8]);
      rd := rd + k*D1;
      rt := rt + k*T1;
      if order >= 2 then
        Rd := 2*g[i,6]*(2*del^2 - g[i,9]*del);
        Rt := 2*g[i,7]*(2*tau^2 - g[i,8]*tau);
        D2 := D1*(D1 - 1) + Rd;
        T2 := T1*(T1 - 1) + Rt;
        rdd := rdd + k*D2;
        rtt := rtt + k*T2;
        rtd := rtd + k*D1*T1;
        if order >= 3 then
          D3 := D1*D2 + (2*D1 - 1)*Rd + 2*g[i,6]*(4*del^2 - g[i,9]*del) - 2*D2;
          T3 := T1*T2 + (2*T1 - 1)*Rt + 2*g[i,7]*(4*tau^2 - g[i,8]*tau) - 2*T2;
          rddd := rddd + k*D3;
          rttt := rttt + k*T3;
          rttd := rttd + k*T2*D1;
          rtdd := rtdd + k*D2*T1;
        end if;
      end if;
    end if;
  end for;

  // divide by delta^m*tau^n
  rd := rd/del;
  rdd := rdd/del^2;
  rddd := rddd/del^3;
  rt := rt/tau;
  rtt := rtt/tau^2;
  rttt := rttt/tau^3;
  rtd := rtd/(del*tau);
  rttd := rttd/(del*tau^2);
  rtdd := rtdd/(del^2*tau);

  // non-analytical terms (Span and Wagner 1996), third derivatives not implemented (as in f_rttt ... f_rddd)
  for i in 1:nNonAna loop
    Psi := exp(-a[i,9]*(delta-1)^2 -a[i,10]*(tau-1)^2);
    Theta := (1-tau) + a[i,8]*((delta-1)^2)^(0.5/a[i,7]);
    Dis := Theta^2 +a[i,11]*((delta-1)^2)^a[i,12];
    Disb := Dis^a[i,6];
    na_r := na_r + a[i,1]*Disb*delta*Psi;
    if order >= 1 then
      Psi_d := -a[i,9]*(2*delta-2)*Psi;
      Theta_d := if (delta-1)<>0 then (a[i,8]*(2*delta-2)*((delta-1)^2)^(0.5/a[i,7])) / (2*a[i,7]*(delta-1)^2) else 0;
      Dis_d := if (delta-1)<>0 then (2*a[i,11]*a[i,12]*((delta-1)^2)^a[i,12]) / (delta-1) +2*Theta*Theta_d else 0;
      Disb_d := if Dis<>0 then a[i,6]*Dis^(a[i,6]-1)*Dis_d else 0;
      Psi_t := -a[i,10]*(2*tau-2)*Psi;
      Dis_t := -2*Theta;
      Disb_t := if Dis<>0 then a[i,6]*Dis^(a[i,6]-1)*Dis_t else 0;
      na_rd := na_rd + a[i,1]*(delta*Disb*Psi_d + delta*Psi*Disb_d + Disb*Psi);
      na_rt := na_rt + a[i,1]*delta*(Disb*Psi_t + Psi*Disb_t);
      if order >= 2 then
        Psi_dd := 2*a[i,9]*(2*a[i,9]*(delta-1)^2-1)*Psi;
        Theta_dd := if (delta-1)<>0 then (a[i,8]*(1/a[i,7]-1)*((delta-1)^2)^(0.5/a[i,7])) / (a[i,7]*(delta-1)^2) else 0;
        Dis_dd := if (delta-1)<>0 then 2*Theta*Theta_dd + 2*Theta_d^2 + (2*((delta-1)^2)^a[i,12]) / (delta-1)^2 *(2*a[i,11]*a[i,12]^2-a[i,11]*a[i,12]) else 0;
        Disb_dd := if (delta-1)<>0 then a[i,6]*(a[i,6]*Dis_d^2 + Dis*Dis_dd - Dis_d^2)*Dis^(a[i,6]-2) else 0;
        Psi_tt := 2*a[i,10]*(2*a[i,10]*(tau-1)^2-1)*Psi;
        Disb_tt := if Dis<>0 then a[i,6]*(a[i,6]*Dis_t^2 + Dis*2 - Dis_t^2)*Dis^(a[i,6]-2)
                   else a[i,6]*(a[i,6]*Dis_t^2 + Dis*2 - Dis_t^2)*Modelica.Constants.inf;
        Psi_td := a[i,9]*a[i,10]*(2*delta-2)*(2*tau-2)*Psi;
        Dis_td := -2*Theta_d;
        Disb_td := if Dis<>0 then a[i,6]/Dis^3*((a[i,6]-1)*Dis^(a[i,6]+1)*Dis_d*Dis_t + Dis^(a[i,6]+2)*Dis_td) else 0;
        na_rdd := na_rdd + a[i,1]*(delta*Disb*Psi_dd + delta*Psi*Disb_dd + 2*delta*Disb_d*Psi_d + 2*Disb*Psi_d + 2*Psi*Disb_d);
        na_rtt := na_rtt + a[i,1]*delta*(Disb*Psi_tt + Psi*Disb_tt + 2*Disb_t*Psi_t);
        na_rtd := na_rtd + a[i,1]*(delta*Disb*Psi_td + delta*Psi*Disb_td + delta*Disb_d*Psi_t + delta*Disb_t*Psi_d + Disb*Psi_t + Psi*Disb_t);
      end if;
    end if;
  end for;
  r := r + na_r;
  rd := rd + na_rd;
  rdd := rdd + na_rdd;
  rt := rt + na_rt;
  rtt := rtt + na_rtt;
  rtd := rtd + na_rtd;

  annotation (Documentation(info="<html>
<p>
Residual part of the dimensionless Helmholtz energy and all its derivatives up to <code>order</code>,
computed in a single loop over all terms. Compared to the separate functions f_r, f_rd, f_rdd, ...,
each term is evaluated only once: the powers delta^d*tau^t and the exponential factors are computed
with one call of exp per term (two for the exponential terms), and all derivatives follow from
multiplying the term value with factors that depend only on d, t, c, eta, beta, gamma and epsilon.
The results agree with the separate functions to rounding precision.
</p>
<p>
The Gaussian terms are implemented with the exponent 2 (columns 4 and 5 of residualGauss), as in the separate functions.
For the non-analytical terms only derivatives up to second order are included, again as in the separate functions.
</p>
</html>"));
end f_r_derivs;
