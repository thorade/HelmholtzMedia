within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function g "returns specifc Gibbs energy g from EoS"
  input HelmholtzDerivs f;
  output SpecificEnergy g;

algorithm
  g := f.T*f.R*((f.i+f.r)+(1+f.delta*f.rd));
end g;
