within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dgdT "returns pressure derivative (dg/dd)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByDensity dgdT;

algorithm
  dgdT := f.T*f.R/f.d*(1+2*f.delta*f.rd + f.delta^2*f.rdd);
end dgdT;
