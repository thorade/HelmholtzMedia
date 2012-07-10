within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_pd
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEnthalpy h "specific enthalpy";
algorithm
  h := specificEnthalpy(setState_pdX(p=p, d=d, phase=phase));
annotation (
  inverse(d=density_ph(p=p, h=h, phase=phase)));
end specificEnthalpy_pd;
