within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_pd "returns specific enthalpy for a given p and d"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
//input ThermodynamicState state;
  output SpecificEnthalpy h "specific enthalpy";

algorithm
  h := Internal.specificEnthalpy_pd_state(p=p, d=d, state=setState_pd(p=p, d=d, phase=phase));

annotation (
  inverse(d=density_ph(p=p, h=h, phase=phase)));
end specificEnthalpy_pd;
