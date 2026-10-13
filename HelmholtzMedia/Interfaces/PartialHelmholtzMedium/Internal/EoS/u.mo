within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.EoS;
function u "returns specifc internal energy u from EoS"
  input HelmholtzDerivs f;
  output SpecificEnergy u;

algorithm
  u := f.T*f.R_s*(f.tau*(f.it + f.rt));
end u;
