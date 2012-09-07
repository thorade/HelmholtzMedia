within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function h "returns specifc enthalpy h from EoS"
  input HelmholtzDerivs f;
  output SpecificEnthalpy h;

algorithm
  h := f.T*f.R*(f.tau*(f.it + f.rt) + (1+f.delta*f.rd));
end h;
