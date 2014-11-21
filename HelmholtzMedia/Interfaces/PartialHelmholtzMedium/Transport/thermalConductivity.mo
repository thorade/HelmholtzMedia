within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity "Return thermal conductivity"
  input ThermodynamicState state;
  output ThermalConductivity lambda;

protected
  ThermalConductivityModel thermalConductivityModel=thermalConductivityCoefficients.thermalConductivityModel;

algorithm
  assert(state.phase <> 2, "thermalConductivity warning: property not defined in two-phase region", level=AssertionLevel.warning);

  if (thermalConductivityModel == ThermalConductivityModel.TC0) then
    // hardcoded models return full thermal conductivity in one equation
    lambda :=thermalConductivity_background(state);
  else
    // composite models
    lambda :=thermalConductivity_dilute(state) + thermalConductivity_background(state) + thermalConductivity_critical(state);
  end if;

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

</html>"));
end thermalConductivity;
