within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Ancillary;
function temperature_pd_Waals
  "start value: temperature from pressure and density, van der Waals equation (gas-like densities)"
  input AbsolutePressure p;
  input Density d;
  output Temperature T;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM "specific gas constant";
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // van der Waals parameters from T_crit and p_crit, as described by Span (2000)
  constant Real a = 27/64*(R_s*R_s*T_crit*T_crit)/p_crit "correction for attraction";
  constant Real b = (R_s*T_crit)/(8*p_crit) "co-volume";
  constant Density d_lim = 0.9/b "upper density limit for the van der Waals equation";

algorithm
  // for d >= d_lim, (1/d - b) approaches zero or becomes negative; the saturation temperature at this (liquid) density
  // is a lower bound of the temperature of compressed liquid with this density
  T := if (d < d_lim) then (p + a*d^2)*(1/d - b)/R_s else Ancillary.saturationTemperature_d(d=d);

  annotation (Documentation(info="<html>
<p>
Temperature from pressure and density, using the van der Waals equation
</p>
<blockquote><pre>T = (p + a&middot;d<sup>2</sup>)&middot;(1/d - b)/R<sub>s</sub></pre></blockquote>
<p>
with a = 27&middot;R<sub>s</sub><sup>2</sup>&middot;T<sub>c</sub><sup>2</sup>/(64&middot;p<sub>c</sub>) and
b = R<sub>s</sub>&middot;T<sub>c</sub>/(8&middot;p<sub>c</sub>), see Span (2000).
This is only meant as a start value for iterations, e.g. in setState_pd.
In the gas and supercritical gas-like region (density up to the critical density) the deviation from the
equation of state is typically 5&nbsp;% to 15&nbsp;%, much smaller than for the ideal gas law.
</p>
<p>
The van der Waals equation is only defined for d &lt; 1/b, which is only about two times the critical density of real fluids.
For d &ge; 0.9/b the saturation temperature at this density (Ancillary.saturationTemperature_d) is returned instead.
</p>
<dl>
<dt>Span, R.</dt>
<dd> <b>Multiparameter Equations of State - An Accurate Source of Thermodynamic Property Data</b>.<br />
     Springer, Berlin, 2000.<br />
     DOI: <a href=\"https://doi.org/10.1007/978-3-662-04092-8\">10.1007/978-3-662-04092-8</a>
</dd>
</dl>
</html>"));
end temperature_pd_Waals;
