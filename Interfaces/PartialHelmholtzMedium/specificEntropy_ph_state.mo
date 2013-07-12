within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEntropy_ph_state
  "returns specific entropy for a given p and h"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific Enthalpy";
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  output SpecificEntropy s "Specific Entropy";

algorithm
  s := specificEntropy(state);

annotation (
  Inline=false,
  LateInline=true);
end specificEntropy_ph_state;
