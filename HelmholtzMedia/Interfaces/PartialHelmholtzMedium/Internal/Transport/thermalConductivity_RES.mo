within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Transport;
function thermalConductivity_RES
  "Returns thermal conductivity (dilute gas and residual contribution, without critical enhancement) from residual entropy scaling"
  input Density d "density";
  input Temperature T "temperature";
  output ThermalConductivity lambda;

protected
  constant Real k_B=Modelica.Constants.k;
  constant Real N_A=Modelica.Constants.N_A;
  constant Real s_plus_min=1e-10 "below this value the residual contribution is neglected (dilute gas)";
  MolarMass MM=fluidConstants[1].molarMass;
  Real[size(thermalConductivityCoefficients.RES_dilute, 1)] c=thermalConductivityCoefficients.RES_dilute;
  Real[4] n=thermalConductivityCoefficients.RES_n;
  Real xi=thermalConductivityCoefficients.RES_xi;
  EoS.HelmholtzDerivs f;
  Real s_plus "reduced residual entropy -s_r/R";
  Real x "scaled reduced residual entropy";
  constant Temperature T_unit=1;
  constant ThermalConductivity lambda_unit=1;
  Real lambda_0 "dilute gas thermal conductivity in W/(m K)";
  Real lambda_plus "plus-scaled residual thermal conductivity";

algorithm
  assert(size(c, 1) > 0, "thermalConductivity_RES: no dilute gas coefficients (RES_dilute) for " + mediumName);
  // dilute gas contribution (Horner scheme, highest power first)
  lambda_0 := 0;
  for i in 1:size(c, 1) loop
    lambda_0 := lambda_0*(T/T_unit) + c[i];
  end for;

  // residual contribution
  f := EoS.setHelmholtzDerivs(d=d, T=T, order=1, includeIdeal=false);
  s_plus := f.r - f.tau*f.rt;
  if (s_plus > s_plus_min) then
    x := s_plus/xi;
    lambda_plus := n[1]*x + n[2]*x^1.5 + n[3]*x^2 + n[4]*x^2.5;
    lambda := lambda_0*lambda_unit + lambda_plus/s_plus^(2/3)*(d/MM*N_A)^(2/3)*k_B*sqrt(k_B*T/(MM/N_A));
  else
    lambda := lambda_0*lambda_unit;
  end if;

  annotation (smoothOrder=2, Documentation(info="<html>
<p>
Thermal conductivity &lambda; = &lambda;<sub>0</sub>(T) + &lambda;<sub>r</sub>, with the residual contribution
</p>
<p>
&lambda;<sub>r</sub> = &lambda;<sub>r</sub><sup>+</sup>/(s<sup>+</sup>)<sup>2/3</sup> &middot; &rho;<sub>N</sub><sup>2/3</sup>&middot;k<sub>B</sub>&middot;(k<sub>B</sub>T/m)<sup>1/2</sup>,
&lambda;<sub>r</sub><sup>+</sup> = n<sub>1</sub>x + n<sub>2</sub>x<sup>1.5</sup> + n<sub>3</sub>x<sup>2</sup> + n<sub>4</sub>x<sup>2.5</sup>,
x = s<sup>+</sup>/&xi;.
</p>
<p>
The critical enhancement (Olchowy and Sengers) is added in
<a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Transport.thermalConductivity\">thermalConductivity</a>.
Coefficients: <code>thermalConductivityCoefficients.RES_dilute</code>, <code>RES_n</code>, <code>RES_xi</code>.
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
<dt>Li, Z., Duan, Y. and Yang, X.</dt>
<dd> <b>Linking Thermal Conductivity to Equations of State Using the Residual Entropy Scaling Theory</b>.<br />
     Ind. Eng. Chem. Res., 63(42):18160-18175, 2024.
</dd>
</dl>
</html>"));
end thermalConductivity_RES;
