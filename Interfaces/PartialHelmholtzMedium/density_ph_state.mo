within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function density_ph_state "returns density for given p and h"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Enthalpy";
  input ThermodynamicState state;
  output Density d "Temperature";

algorithm
  d := density(state);

annotation (
  derivative(noDerivative=state)=density_ph_der,
  Inline=false,
  LateInline=true);
end density_ph_state;
