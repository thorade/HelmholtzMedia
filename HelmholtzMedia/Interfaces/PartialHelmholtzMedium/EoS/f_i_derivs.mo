within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_i_derivs
  "ideal part of the dimensionless Helmholtz energy and its tau-derivatives up to the given order, in one pass over all terms"
  input Real delta "reduced density";
  input Real tau "inverse reduced temperature";
  input Integer order(min=0, max=3)=2 "highest order of the derivatives to compute";
  output Real i "f_i";
  output Real it "(df_i/dtau)@delta=const";
  output Real itt "(d2f_i/dtau2)@delta=const";
  output Real ittt "(d3f_i/dtau3)@delta=const";

protected
  constant Integer nLog=size(helmholtzCoefficients.idealLog, 1);
  constant Integer nPower=size(helmholtzCoefficients.idealPower, 1);
  constant Integer nEinstein=size(helmholtzCoefficients.idealEinstein, 1);
  constant Integer nCosh=size(helmholtzCoefficients.idealCosh, 1);
  constant Integer nSinh=size(helmholtzCoefficients.idealSinh, 1);

  constant Real[nLog, 2] l=helmholtzCoefficients.idealLog;
  constant Real[nPower, 2] p=helmholtzCoefficients.idealPower;
  constant Real[nEinstein, 2] e=helmholtzCoefficients.idealEinstein;
  constant Real[nCosh, 2] c=helmholtzCoefficients.idealCosh;
  constant Real[nSinh, 2] s=helmholtzCoefficients.idealSinh;

  Real lnt "log(tau)";
  Real w "value of the current term";
  Real E "exp(b*tau) (Einstein terms)";
  Real th "tanh(b*tau) (cosh and sinh terms)";

algorithm
  i := 0;
  it := 0;
  itt := 0;
  ittt := 0;
  if (delta > 0) and (tau > 0) then
    lnt := log(tau);
    i := log(delta);
    for j in 1:nLog loop
      // l1*log(tau^l2)
      i := i + l[j,1]*l[j,2]*lnt;
      if order >= 1 then
        it := it + l[j,1]*l[j,2]/tau;
        if order >= 2 then
          itt := itt - l[j,1]*l[j,2]/tau^2;
          if order >= 3 then
            ittt := ittt + 2*l[j,1]*l[j,2]/tau^3;
          end if;
        end if;
      end if;
    end for;
    for j in 1:nPower loop
      // p1*tau^p2
      w := p[j,1]*exp(p[j,2]*lnt);
      i := i + w;
      if order >= 1 then
        it := it + w*p[j,2]/tau;
        if order >= 2 then
          itt := itt + w*p[j,2]*(p[j,2] - 1)/tau^2;
          if order >= 3 then
            ittt := ittt + w*p[j,2]*(p[j,2] - 1)*(p[j,2] - 2)/tau^3;
          end if;
        end if;
      end if;
    end for;
    for j in 1:nEinstein loop
      // e1*log(1 - exp(e2*tau))
      E := exp(e[j,2]*tau);
      i := i + e[j,1]*log(1 - E);
      if order >= 1 then
        it := it - e[j,1]*e[j,2]*E/(1 - E);
        if order >= 2 then
          itt := itt - e[j,1]*e[j,2]^2*E/(1 - E)^2;
          if order >= 3 then
            ittt := ittt + e[j,1]*e[j,2]^3*E*(E + 1)/(E - 1)^3;
          end if;
        end if;
      end if;
    end for;
    for j in 1:nCosh loop
      // -c1*log(|cosh(c2*tau)|)
      i := i - c[j,1]*log(abs(cosh(c[j,2]*tau)));
      if order >= 1 then
        th := tanh(c[j,2]*tau);
        it := it - c[j,1]*c[j,2]*th;
        if order >= 2 then
          itt := itt - c[j,1]*c[j,2]^2*(1 - th^2);
          if order >= 3 then
            ittt := ittt - 2*c[j,1]*c[j,2]^3*th*(th^2 - 1);
          end if;
        end if;
      end if;
    end for;
    for j in 1:nSinh loop
      // s1*log(|sinh(s2*tau)|)
      i := i + s[j,1]*log(abs(sinh(s[j,2]*tau)));
      if order >= 1 then
        th := tanh(s[j,2]*tau);
        it := it + s[j,1]*s[j,2]/th;
        if order >= 2 then
          itt := itt - s[j,1]*s[j,2]^2*(1/th^2 - 1);
          if order >= 3 then
            ittt := ittt + 2*s[j,1]*s[j,2]^3/th*(1/th^2 - 1);
          end if;
        end if;
      end if;
    end for;
  else
    i := -Modelica.Constants.inf;
  end if;

  annotation (Documentation(info="<html>
<p>
Ideal part of the dimensionless Helmholtz energy and its tau-derivatives up to <code>order</code>,
computed in a single loop over all terms (log(tau) and the exponential or hyperbolic functions once per term).
The results agree with the separate functions f_i, f_it, f_itt and f_ittt to rounding precision.
</p>
</html>"));
end f_i_derivs;
