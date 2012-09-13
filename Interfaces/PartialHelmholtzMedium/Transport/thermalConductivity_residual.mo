within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity_residual
  "Return thermal conductivity residual contribution"
  input ThermodynamicState state;
  output ThermalConductivity lambda_r;

protected
  MolarMass MM = fluidConstants[1].molarMass;
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_residual=fluidConstants[1].molarMass/
      thermalConductivityCoefficients.reducingMolarVolume_residual;
  Real delta "reduced density";

  Temperature T_red_residual=thermalConductivityCoefficients.reducingTemperature_residual;
  Real tau "reduced temperature";

  // coeffs for residual contribution
  Integer nResidual = size(thermalConductivityCoefficients.lambda_r_coeffs, 1);
  Real[nResidual, 4] B= thermalConductivityCoefficients.lambda_r_coeffs;

  Real lambda_red_residual=thermalConductivityCoefficients.reducingThermalConductivity_residual;

algorithm
  // residual contribution; RefProp uses the name background contribution
  tau := state.T/T_red_residual;
  delta := state.d/d_red_residual;
  lambda_r := sum((B[i, 1]*tau^B[i, 2])*(delta)^B[i, 3] for i in 1:nResidual);
  lambda_r := lambda_r*lambda_red_residual;

  /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        d = " + String(state.d) + " and T = " + String(state.T));
  Modelica.Utilities.Streams.print(" lambda_r = " + String(lambda_r));
  Modelica.Utilities.Streams.print("===========================================");
  */

  annotation (Documentation(info="<html>
  <p>
The thermal conductivity (TC) is split into three parts: ideal gas TC lamda_0, residual TC lambda_r and critical TC enhancement lambda_c.
Sometimes the residual TC is again split into two parts.
This allows to develop functions for each contribution seperately.
The sum of ideal gas TC and residual TC is called background TC.
Ideal gas TC depends on Temperature only and can be modelled by a quadratic function.
Residual TC is also modeled by a polynominal.
At the critical point TC becomes infinite; TC is enhanced for a large region around the critical point.
The critical enhancement can be described by various alternative approaches.
Here, the simplified approach as suggested by Olchowy and Sengers is implemented.

Special thanks go to Eric W. Lemmon for answering all my emails 
and programming a special version of RefProp that outputs also intermediate values.

</p>
<dl>
<dt>Olchowy, G.A. and Sengers, J.V</dt>
<dd> <b>A simplified representation for the thermal conductivity of fluids in the critical region</b>.<br>
     International Journal of Thermophysics (1998) 10, 417-426.<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
end thermalConductivity_residual;
