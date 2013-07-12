within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEntropy_dT "return specific enthalpy for given d and T"
  extends Modelica.Icons.Function;
  input Density d "Density";
  input Temperature T "Temperature";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
//input ThermodynamicState state;
  output SpecificEntropy s "specific entropy";

algorithm
  s := specificEntropy(setState_dTX(d=d,T=T,phase=phase));

annotation (
  inverse(d=density_Ts(T=T, s=s, phase=phase)));
end specificEntropy_dT;
