within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_Ts "returns pressure for given T and s"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Specific entropy";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
//input ThermodynamicState state;
  output AbsolutePressure p "Pressure";

algorithm
  p := pressure(setState_Ts(T=T, s=s, phase=phase));

annotation (
  inverse(s=specificEntropy_pT(p=p, T=T, phase=phase),
          T=temperature_ps(p=p, s=s, phase=phase)));
end pressure_Ts;
