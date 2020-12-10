within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function p "returns pressure p from EoS"
  input HelmholtzDerivs f;
  output AbsolutePressure p;

algorithm
  p := f.d*f.T*f.R_s*(1+f.delta*f.rd);
end p;
