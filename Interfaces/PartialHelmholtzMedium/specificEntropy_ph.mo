within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEntropy_ph "returns specific entropy for a given p and h"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific Enthalpy";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEntropy s "Specific Entropy";

algorithm
  s := specificEntropy(setState_phX(p=p, h=h, phase=1));

annotation (
  inverse(h=specificEnthalpy_ps(p=p, s=s, phase=phase)));
end specificEntropy_ph;
