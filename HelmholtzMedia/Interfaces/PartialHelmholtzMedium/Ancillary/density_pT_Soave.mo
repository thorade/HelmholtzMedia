within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function density_pT_Soave "solve Redlich-Kwong-Soave for density"
  input Temperature T;
  input AbsolutePressure p;
  input AbsolutePressure psat=Ancillary.saturationPressure_T(T=T)
    "selects the liquid (p>psat) or vapour root if there are several physical roots";
  output Density d;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // RKS: Redlich-Kwong-Soave (see Soave 1972)
  constant Real omega = fluidConstants[1].acentricFactor;
  constant Real m = 0.480 + 1.574*omega - 0.176*omega^2;
  Real a = 0.42747*R_s^2*T_crit^2/p_crit*(1 + m*(1 - sqrt(T/T_crit)))^2;
  Real b = 0.08664*R_s*T_crit/p_crit;
  Real A = a*p/(R_s^2*T^2);
  Real B = b*p/(R_s*T);

  // cubic in the compressibility factor: Z^3 + a2*Z^2 + a1*Z + a0 = 0;
  // the constant term is proportional to p^2, the polynomial stays well-conditioned at low pressure
  Real a2 = -1;
  Real a1 = A - B - B^2;
  Real a0 = -A*B;
  Real x_infl = -a2/3 "inflection point";
  Real D = a2^2 - 3*a1 "the cubic has a maximum and a minimum if D>0";
  Real x "first root";
  Real fx;
  Real dfx;
  Real d2fx;
  Real dx = 1;
  Real c1 "deflated quadratic: Z^2 + c1*Z + c0";
  Real c0;
  Real disc;
  Real q;
  Real Z[3] "real roots";
  Integer nRoots = 1;
  Real Z_sel = 0 "selected physical root";
  Integer iter = 0;
  constant Integer iter_max = 100;

algorithm
  // all real roots following Deiters and Macias-Salinas (2014): first root by Halley's method,
  // started at the inflection point or at a Laguerre-Nair-Samuelson bound, so that no extremum lies
  // between start value and root; then deflation and the analytical solution of the quadratic
  if (D > 0) then
    x := if (((x_infl + a2)*x_infl + a1)*x_infl + a0 > 0) then x_infl - 2/3*sqrt(D) else x_infl + 2/3*sqrt(D);
  else
    x := x_infl;
  end if;
  while (iter < iter_max) and (abs(dx) > 1e-15*max(1, abs(x))) loop
    iter := iter + 1;
    fx := ((x + a2)*x + a1)*x + a0;
    dfx := (3*x + 2*a2)*x + a1;
    d2fx := 6*x + 2*a2;
    dx := if (abs(dfx*dfx - 0.5*fx*d2fx) > 0) then fx*dfx/(dfx*dfx - 0.5*fx*d2fx) else 0;
    x := x - dx;
  end while;

  // deflation: if the first root is the largest one, use Vieta's rules, which avoid the cancellation
  // in c0 = c1*x + a1 that destroys the small (liquid) roots at low pressure
  if (abs(x)^3 >= abs(a0)) then
    c0 := -a0/x;
    c1 := -(a1 - c0)/x;
  else
    c1 := x + a2;
    c0 := c1*x + a1;
  end if;
  Z := {x, 0, 0};
  disc := c1^2/4 - c0;
  if (disc >= 0) then
    // quadratic formula without cancellation
    q := -(c1/2 + (if c1 >= 0 then 1 else -1)*sqrt(disc));
    if (abs(q) > 0) then
      Z[2] := q;
      Z[3] := c0/q;
      nRoots := 3;
    end if;
  end if;

  // physical roots have a specific volume larger than the covolume b, i.e. Z>B;
  // of several physical roots, the liquid uses the smallest and the vapour the largest
  for i in 1:nRoots loop
    if (Z[i] > B) then
      if (Z_sel <= 0) then
        Z_sel := Z[i];
      elseif (p > psat) then
        Z_sel := min(Z_sel, Z[i]);
      else
        Z_sel := max(Z_sel, Z[i]);
      end if;
    end if;
  end for;
  if (Z_sel <= 0) then
    Z_sel := x;
  end if;
  d := p/(R_s*T*Z_sel);

  annotation (Documentation(info="<html>
<p>Density from the Redlich-Kwong-Soave equation of state, used as start value for the iterative solution of the
Helmholtz EoS (see Span 2000, section 3.3.1). All real roots of the cubic equation in the compressibility factor Z
are determined with the algorithm of Deiters and Mac&iacute;as-Salinas (2014); of those, only roots with Z&gt;B
(specific volume larger than the covolume b) are physical. If there are several physical roots, the liquid root
(smallest Z) is used for p&gt;psat and the vapour root (largest Z) otherwise.</p>
<dl>
<dt>Soave, G.</dt>
<dd><b>Equilibrium constants from a modified Redlich-Kwong equation of state</b><br />
     Chemical Engineering Science 27, 1197-1203 (1972)</dd>
<dt>Deiters, U. K. and Mac&iacute;as-Salinas, R.</dt>
<dd><b>Calculation of Densities from Cubic Equations of State: Revisited</b><br />
     Ind. Eng. Chem. Res. 53, 2529-2536 (2014)<br />
     DOI: <a href=\"http://dx.doi.org/10.1021/ie4038664\">10.1021/ie4038664</a></dd>
</dl>
</html>"));
end density_pT_Soave;
