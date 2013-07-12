within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function temperature_pd "returns temperature for given p and d"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
//input ThermodynamicState state;
  output Temperature T "Temperature";

algorithm
  T := temperature(setState_pd(p=p, d=d, phase=phase));

annotation (
  inverse(p=pressure_dT(d=d, T=T, phase=phase),
          d=density_pT(p=p, T=T, phase=phase)));
end temperature_pd;
