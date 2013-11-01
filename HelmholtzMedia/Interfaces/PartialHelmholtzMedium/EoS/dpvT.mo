within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dpvT "returns pressure derivative (dp/dv)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerPressureByVolume dpvT;

algorithm
  dpvT := -f.d*f.d*f.T*f.R*(1 + 2*f.delta*f.rd + f.delta*f.delta*f.rdd);
end dpvT;
