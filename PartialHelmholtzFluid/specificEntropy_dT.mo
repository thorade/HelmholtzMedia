within HelmholtzFluids.PartialHelmholtzFluid;
function specificEntropy_dT "Return specific enthalpy from d and T"
  extends Modelica.Icons.Function;

  input Density d "Density";
  input Temperature T "Temperature";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEntropy s "specific entropy";

algorithm
  s := specificEntropy(setState_dTX(d=d,T=T,X=fill(0, 0),phase=phase));
end specificEntropy_dT;
