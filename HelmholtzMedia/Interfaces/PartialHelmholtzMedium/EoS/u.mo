within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function u "returns specifc internal energy u from EoS"
  input HelmholtzDerivs f;
  output SpecificEnergy u;

algorithm
  u := f.T*f.R*(f.tau*(f.it + f.rt));
end u;
