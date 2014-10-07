within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function surfaceTension "Return surface tension sigma in the two phase region"
    input SaturationProperties sat;
    output SurfaceTension sigma;

protected
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real[size(surfaceTensionCoefficients.coeffs, 1)] a=surfaceTensionCoefficients.coeffs[:, 1];
  Real[size(surfaceTensionCoefficients.coeffs, 1)] n=surfaceTensionCoefficients.coeffs[:, 2];
  Real X "reduced temperature difference";

algorithm
  X := (T_crit - sat.Tsat)/T_crit;
  sigma := sum(a[i]*X^n[i] for i in 1:size(a, 1));

  annotation (Documentation(info="<html>
  <p>
This is an implementation of the model as suggested by Somayajulu, G.R.,
which is an extension of the van der Waals surface tension correlation. 
The extended version has up to three terms with two parameters each.
This algorithm uses T only, some mixture models also use liquid and vapour density.
</p>
<dl>
<dt>Somayajulu, G.R.</dt>
<dd> <b>A generalized equation for surface tension from the triple point to the critical point</b>.<br>
     International Journal of Thermophysics (1988) 9, 559-566.<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF00503154\">10.1007/BF00503154</a>
</dd>
<dt>Van der Waals, J.D.</dt>
<dd> <b>Thermodynamische Theorie der Kapillarit&auml;t unter Voraussetzung stetiger Dichte&auml;nderung</b>.<br>
     Zeitschrift f&uuml;r Physikalische Chemie (1894) 13, 657-725.
</dd>
</dl>
</html>"));
end surfaceTension;
