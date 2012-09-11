within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dpdT "returns pressure derivative (dp/dd)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerPressureByDensity dpdT;

algorithm
  dpdT := f.T*f.R*(1 + 2*f.delta*f.rd + f.delta*f.delta*f.rdd);
end dpdT;
