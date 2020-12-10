within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dgdT "returns Gibbs energy derivative (dg/dd)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByDensity dgdT;

algorithm
  dgdT := f.T*f.R_s/f.d*(1+2*f.delta*f.rd + f.delta*f.delta*f.rdd);
end dgdT;
