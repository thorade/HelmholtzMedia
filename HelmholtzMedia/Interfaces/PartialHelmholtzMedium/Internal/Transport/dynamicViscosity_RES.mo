within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Transport;
function dynamicViscosity_RES
  "Returns dynamic viscosity (in micro Pa s) from residual entropy scaling"
  input Density d "density";
  input Temperature T "temperature";
  output Real eta "dynamic viscosity in micro Pa s";

protected
  constant Real k_B=Modelica.Constants.k;
  constant Real N_A=Modelica.Constants.N_A;
  constant Real micro=1e-6;
  constant Temperature T_unit=1;
  constant Real s_plus_min=1e-10 "below this value the residual contribution is neglected (dilute gas)";
  MolarMass MM=fluidConstants[1].molarMass;
  Real[size(dynamicViscosityCoefficients.RES_dilute, 1)] c=dynamicViscosityCoefficients.RES_dilute;
  Real[3] n=dynamicViscosityCoefficients.RES_n;
  Real xi=dynamicViscosityCoefficients.RES_xi;
  EoS.HelmholtzDerivs f;
  Real s_plus "reduced residual entropy -s_r/R";
  Real x "scaled reduced residual entropy";
  Real eta_0 "dilute gas viscosity in micro Pa s";
  Real eta_plus "plus-scaled residual viscosity";

algorithm
  assert(size(c, 1) > 0, "dynamicViscosity_RES: no dilute gas coefficients (RES_dilute) for " + mediumName);
  // dilute gas contribution (Horner scheme, highest power first)
  eta_0 := 0;
  for i in 1:size(c, 1) loop
    eta_0 := eta_0*(T/T_unit) + c[i];
  end for;

  // residual contribution
  f := EoS.setHelmholtzDerivs(d=d, T=T, order=1, includeIdeal=false);
  s_plus := f.r - f.tau*f.rt;
  if (s_plus > s_plus_min) then
    x := s_plus/xi;
    eta_plus := exp(n[1]*x^1.8 + n[2]*x^2.4 + n[3]*x^2.8) - 1;
    eta := eta_0 + eta_plus/s_plus^(2/3)*(d/MM*N_A)^(2/3)*sqrt(MM/N_A*k_B*T)/micro;
  else
    eta := eta_0;
  end if;

  annotation (smoothOrder=2, Documentation(info="<html>
<p>
Dynamic viscosity &eta; = &eta;<sub>0</sub>(T) + &eta;<sub>r</sub>, with the residual contribution
</p>
<p>
&eta;<sub>r</sub> = &eta;<sub>r</sub><sup>+</sup>/(s<sup>+</sup>)<sup>2/3</sup> &middot; &rho;<sub>N</sub><sup>2/3</sup>&middot;(m&middot;k<sub>B</sub>&middot;T)<sup>1/2</sup>,
ln(&eta;<sub>r</sub><sup>+</sup> + 1) = n<sub>1</sub>x<sup>1.8</sup> + n<sub>2</sub>x<sup>2.4</sup> + n<sub>3</sub>x<sup>2.8</sup>,
x = s<sup>+</sup>/&xi;,
</p>
<p>
where &rho;<sub>N</sub> is the number density and m the mass of one molecule.
Coefficients: <code>dynamicViscosityCoefficients.RES_dilute</code>, <code>RES_n</code>, <code>RES_xi</code>.
</p>
<p>
Residual entropy scaling (RES): the residual part of the transport property is a function of the reduced residual entropy
s<sup>+</sup> = -s<sub>r</sub>/R = &alpha;<sup>r</sup> - &tau;&middot;(&part;&alpha;<sup>r</sup>/&part;&tau;)<sub>&delta;</sub>,
which is computed from the Helmholtz equation of state of the fluid.
The dilute gas contribution is a polynomial in T. Below s<sup>+</sup> = 1e-10 (dilute gas) the residual contribution is neglected.
</p>
<p>
The RES parameters were fitted with the residual entropy of the equations of state in REFPROP 10.
If the fluid uses a different equation of state, the results deviate accordingly.
</p>
<dl>
<dt>Martinek, V., Bell, I.H., Herzog, R., Richter, M. and Yang, X.</dt>
<dd> <b>Entropy Scaling of Viscosity IV - Application to 124 Industrially Important Fluids</b>.<br />
     J. Chem. Eng. Data, 2025.
</dd>
</dl>
</html>"));
end dynamicViscosity_RES;
