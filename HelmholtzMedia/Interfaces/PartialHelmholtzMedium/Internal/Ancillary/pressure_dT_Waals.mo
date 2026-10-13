within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Ancillary;
function pressure_dT_Waals
  "start value: pressure from density and temperature, van der Waals equation (gas-like densities)"
  input Density d;
  input Temperature T;
  output AbsolutePressure p;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM "specific gas constant";
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // van der Waals parameters from T_crit and p_crit, as described by Span (2000)
  constant Real a = 27/64*(R_s*R_s*T_crit*T_crit)/p_crit "correction for attraction";
  constant Real b = (R_s*T_crit)/(8*p_crit) "co-volume";
  constant Density d_lim = 0.9/b "upper density limit for the van der Waals equation";

  AbsolutePressure p_lim "pressure at d_lim";
  Real dpdd_lim "(dp/dd)@T at d_lim";

algorithm
  if (d < d_lim) then
    p := R_s*T/(1/d - b) - a*d^2;
  else
    // the van der Waals pressure diverges at d=1/b; continue linearly beyond d_lim instead
    p_lim := R_s*T/(1/d_lim - b) - a*d_lim^2;
    dpdd_lim := R_s*T/(1 - b*d_lim)^2 - 2*a*d_lim;
    p := p_lim + dpdd_lim*(d - d_lim);
  end if;

  annotation (Documentation(info="<html>
<p>
Pressure from density and temperature, using the van der Waals equation
</p>
<blockquote><pre>p = R<sub>s</sub>&middot;T/(1/d - b) - a&middot;d<sup>2</sup></pre></blockquote>
<p>
with a = 27&middot;R<sub>s</sub><sup>2</sup>&middot;T<sub>c</sub><sup>2</sup>/(64&middot;p<sub>c</sub>) and
b = R<sub>s</sub>&middot;T<sub>c</sub>/(8&middot;p<sub>c</sub>), see Span (2000).
This is only meant as a rough estimate or start value, for gas-like densities.
</p>
<p>
The van der Waals equation is only defined for d &lt; 1/b, which is only about two times the critical density of real fluids.
For d &ge; 0.9/b the pressure is extrapolated linearly from d = 0.9/b, which keeps it finite and continuous,
but it is not a good estimate for compressed liquid.
</p>
<dl>
<dt>Span, R.</dt>
<dd> <b>Multiparameter Equations of State - An Accurate Source of Thermodynamic Property Data</b>.<br />
     Springer, Berlin, 2000.<br />
     DOI: <a href=\"https://doi.org/10.1007/978-3-662-04092-8\">10.1007/978-3-662-04092-8</a>
</dd>
</dl>
</html>"));
end pressure_dT_Waals;
