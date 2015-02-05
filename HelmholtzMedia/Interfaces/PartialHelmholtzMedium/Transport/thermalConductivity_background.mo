within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity_background
  "Return thermal conductivity background/residual  contribution"
  input ThermodynamicState state;
  output ThermalConductivity lambda_b;

protected
  ThermalConductivityModel thermalConductivityModel=thermalConductivityCoefficients.thermalConductivityModel;
  MolarMass MM = fluidConstants[1].molarMass;
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_background=fluidConstants[1].molarMass/thermalConductivityCoefficients.reducingMolarVolume_background;
  Real delta "reduced density";
  Temperature T_red_background=thermalConductivityCoefficients.reducingTemperature_background;
  Real tau "reduced temperature";

  // coeffs for background contribution
  Integer nBackground = size(thermalConductivityCoefficients.lambda_b_coeffs, 1);
  Real[nBackground, 4] B= thermalConductivityCoefficients.lambda_b_coeffs;

  Real lambda_red_background=thermalConductivityCoefficients.reducingThermalConductivity_background;

algorithm
  if (thermalConductivityModel == ThermalConductivityModel.TC0) then
    // hardcoded models, use mediumName to distinguish further
    if mediumName == "helium" then
    end if;
  elseif (thermalConductivityModel == ThermalConductivityModel.TC1) then
    tau := state.T/T_red_background;
    delta := state.d/d_red_background;
    lambda_b := sum((B[i, 1]*tau^B[i, 2])*(delta)^B[i, 3] for i in 1:nBackground);
    lambda_b := lambda_b*lambda_red_background;
  elseif (thermalConductivityModel == ThermalConductivityModel.TC2) then
    assert(false, "ThermalconductivityModel TC2 not yet implemented");
  else
    assert(false, "unknown ThermalconductivityModel");
  end if;

  /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        d = " + String(state.d) + " and T = " + String(state.T));
  Modelica.Utilities.Streams.print(" lambda_b = " + String(lambda_b));
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
<dd> <b>A simplified representation for the thermal conductivity of fluids in the critical region</b>.<br />
     International Journal of Thermophysics (1998) 10, 417-426.<br />
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
end thermalConductivity_background;
