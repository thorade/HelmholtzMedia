within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function temperature_ph_state "returns temperature for given p and h"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Enthalpy";
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  output Temperature T "Temperature";

algorithm
  T := temperature(state);

annotation (
  Inline=false,
  LateInline=true,
  inverse(h=specificEnthalpy_pT_state(p=p, T=T, state=state)),
  derivative(noDerivative=state)=temperature_ph_der);
end temperature_ph_state;
