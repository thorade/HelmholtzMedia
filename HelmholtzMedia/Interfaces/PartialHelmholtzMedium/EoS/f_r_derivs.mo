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
  constant Real gaussExponentDeviation = sum(abs(g[i,4] - 2) + abs(g[i,5] - 2) for i in 1:nGauss)
    "zero if all Gaussian terms use the exponent 2 (evaluated at translation, no cost at run time)";

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

  // non-analytical terms
  Real u "delta-1";
  Real v "tau-1";
  Real e "exponent of |delta-1|";
  Real h "|delta-1|^(1/beta)";
  Real h1 "dh/ddelta";
  Real h2 "d2h/ddelta2";
  Real h3 "d3h/ddelta3";
  Real q "|delta-1|^(2*a)";
  Real q1 "dq/ddelta";
  Real q2 "d2q/ddelta2";
  Real q3 "d3q/ddelta3";
  Real Theta;
  Real Theta_d;
  Real Theta_dd;
  Real Dis "distance function Delta";
  Real Dis_d;
  Real Dis_t;
  Real Dis_dd;
  Real Dis_td;
  Real Dis_ddd;
  Real Dis_tdd;
  Real Psi;
  Real F "Delta^b";
  Real f1;
  Real f2;
  Real f3;
  Real F_d;
  Real F_t;
  Real F_dd;
  Real F_tt;
  Real F_td;
  Real F_ddd;
  Real F_ttt;
  Real F_tdd;
  Real F_ttd;
  Real p1;
  Real p2;
  Real p3;
  Real s1;
  Real s2;
  Real s3;
  Real G_d "derivatives of G = F*Psi, divided by Psi";
  Real G_t;
  Real G_dd;
  Real G_tt;
  Real G_td;
  Real na_r=0;
  Real na_rd=0;
  Real na_rdd=0;
  Real na_rt=0;
  Real na_rtt=0;
  Real na_rtd=0;
  Real na_rttt=0;
  Real na_rttd=0;
  Real na_rtdd=0;
  Real na_rddd=0;

algorithm
  assert(gaussExponentDeviation == 0,
    "f_r_derivs: only Gaussian terms with exponent 2 (columns 4 and 5 of residualGauss) are implemented");
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

  // non-analytical terms n*Delta^b*delta*Psi (Span and Wagner 1996), with
  //   Delta = Theta^2 + B*((delta-1)^2)^a,  Theta = (1-tau) + A*((delta-1)^2)^(1/(2*beta)),  Psi = exp(-C*(delta-1)^2 - D*(tau-1)^2).
  // The derivatives of the building blocks are combined with the chain rule (F = Delta^b) and the Leibniz rule (G = F*Psi).
  // In the derivatives of |delta-1|^e the exponents e = 1/beta and e = 2*a are larger than 3 for all known coefficient sets,
  // so these derivatives vanish at delta = 1. At the critical point (Delta = 0) the derivatives of Delta^b diverge; they are set to zero there.
  u := delta - 1;
  v := tau - 1;
  for i in 1:nNonAna loop
    // h = |u|^(1/beta) and q = |u|^(2*a), with derivatives w.r.t. delta
    if u <> 0 then
      e := 1/a[i,7];
      h := abs(u)^e;
      h1 := e*h/u;
      h2 := (e - 1)*h1/u;
      h3 := (e - 2)*h2/u;
      e := 2*a[i,12];
      q := abs(u)^e;
      q1 := e*q/u;
      q2 := (e - 1)*q1/u;
      q3 := (e - 2)*q2/u;
    else
      h := 0;
      h1 := 0;
      h2 := 0;
      h3 := 0;
      q := 0;
      q1 := 0;
      q2 := 0;
      q3 := 0;
    end if;
    Theta := -v + a[i,8]*h;
    Dis := Theta^2 + a[i,11]*q;
    F := if Dis > 0 then Dis^a[i,6] else 0;
    Psi := exp(-a[i,9]*u^2 - a[i,10]*v^2);
    k := a[i,1]*Psi;
    na_r := na_r + k*delta*F;
    if order >= 1 then
      Theta_d := a[i,8]*h1;
      Dis_d := 2*Theta*Theta_d + a[i,11]*q1;
      Dis_t := -2*Theta;
      f1 := if Dis > 0 then a[i,6]*F/Dis else 0 "dF/dDelta";
      p1 := -2*a[i,9]*u "Psi_delta/Psi";
      s1 := -2*a[i,10]*v "Psi_tau/Psi";
      F_d := f1*Dis_d;
      F_t := f1*Dis_t;
      G_d := F_d + F*p1;
      G_t := F_t + F*s1;
      na_rd := na_rd + k*(F + delta*G_d);
      na_rt := na_rt + k*delta*G_t;
      if order >= 2 then
        Theta_dd := a[i,8]*h2;
        Dis_dd := 2*Theta_d^2 + 2*Theta*Theta_dd + a[i,11]*q2;
        Dis_td := -2*Theta_d;
        // Dis_tt = 2
        f2 := if Dis > 0 then (a[i,6] - 1)*f1/Dis else 0 "d2F/dDelta2";
        p2 := p1^2 - 2*a[i,9];
        s2 := s1^2 - 2*a[i,10];
        F_dd := f2*Dis_d^2 + f1*Dis_dd;
        F_tt := f2*Dis_t^2 + 2*f1;
        F_td := f2*Dis_d*Dis_t + f1*Dis_td;
        G_dd := F_dd + 2*F_d*p1 + F*p2;
        G_tt := F_tt + 2*F_t*s1 + F*s2;
        G_td := F_td + F_d*s1 + F_t*p1 + F*p1*s1;
        na_rdd := na_rdd + k*(2*G_d + delta*G_dd);
        na_rtt := na_rtt + k*delta*G_tt;
        na_rtd := na_rtd + k*(G_t + delta*G_td);
        if order >= 3 then
          Dis_ddd := 6*Theta_d*Theta_dd + 2*Theta*a[i,8]*h3 + a[i,11]*q3;
          Dis_tdd := -2*Theta_dd;
          // Dis_ttt = Dis_ttd = 0
          f3 := if Dis > 0 then (a[i,6] - 2)*f2/Dis else 0 "d3F/dDelta3";
          p3 := p1^3 - 6*a[i,9]*p1;
          s3 := s1^3 - 6*a[i,10]*s1;
          F_ddd := f3*Dis_d^3 + 3*f2*Dis_d*Dis_dd + f1*Dis_ddd;
          F_ttt := f3*Dis_t^3 + 6*f2*Dis_t;
          F_tdd := f3*Dis_d^2*Dis_t + f2*(Dis_dd*Dis_t + 2*Dis_td*Dis_d) + f1*Dis_tdd;
          F_ttd := f3*Dis_d*Dis_t^2 + f2*(2*Dis_d + 2*Dis_td*Dis_t);
          na_rddd := na_rddd + k*(3*G_dd + delta*(F_ddd + 3*F_dd*p1 + 3*F_d*p2 + F*p3));
          na_rttt := na_rttt + k*delta*(F_ttt + 3*F_tt*s1 + 3*F_t*s2 + F*s3);
          na_rttd := na_rttd + k*(G_tt + delta*(F_ttd + F_tt*p1 + 2*(F_td + F_t*p1)*s1 + (F_d + F*p1)*s2));
          na_rtdd := na_rtdd + k*(2*G_td + delta*(F_tdd + F_dd*s1 + 2*(F_td + F_d*s1)*p1 + (F_t + F*s1)*p2));
        end if;
      end if;
    end if;
  end for;
  r := r + na_r;
  rd := rd + na_rd;
  rdd := rdd + na_rdd;
  rt := rt + na_rt;
  rtt := rtt + na_rtt;
  rtd := rtd + na_rtd;
  rttt := rttt + na_rttt;
  rttd := rttd + na_rttd;
  rtdd := rtdd + na_rtdd;
  rddd := rddd + na_rddd;

  annotation (Documentation(info="<html>
<p>
Residual part of the dimensionless Helmholtz energy and all its derivatives up to <code>order</code>,
computed in a single loop over all terms. Each term is evaluated only once: the powers delta^d*tau^t and the
exponential factors are computed with one call of exp per term (two for the exponential terms), and all derivatives
follow from multiplying the term value with factors D_m and T_n that depend only on d, t, c, eta, beta, gamma and epsilon.
This factored form also avoids the cancellation of fully expanded derivative expressions.
</p>
<p>
The Gaussian terms are implemented with the exponent 2 only (columns 4 and 5 of residualGauss); an assert checks this.
The non-analytical terms (Span and Wagner 1996, used for carbon dioxide) are included up to third order.
</p>
<p>
The factors of the polynomial, exponential and Gaussian terms were verified symbolically (sympy);
the non-analytical terms were verified against high-precision numerical differentiation (mpmath, 60 digits) and
against <a href=\"https://github.com/usnistgov/teqp\">teqp</a> (automatic differentiation),
see Examples.Validation.HelmholtzDerivsNonAnalytical.
</p>
</html>"));
end f_r_derivs;
