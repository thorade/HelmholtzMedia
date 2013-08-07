within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function s "returns specifc entropy s from EoS"
  input HelmholtzDerivs f;
  output SpecificEntropy s;

algorithm
  s := f.R*(f.tau*(f.it + f.rt) - f.i - f.r);
end s;
