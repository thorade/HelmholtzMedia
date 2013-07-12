within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function density_Ts "returns density for given T and s"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Specific entropy";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
//input ThermodynamicState state;
  output Density d "Density";

algorithm
  d := density(setState_Ts(T=T, s=s, phase=phase));

annotation (
  inverse=specificEntropy_dT(d=d, T=T, phase=phase));
end density_Ts;
