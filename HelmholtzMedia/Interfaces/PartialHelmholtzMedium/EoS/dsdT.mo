within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dsdT "returns entropy derivative (ds/dd)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerEntropyByDensity dsdT;

algorithm
  dsdT := -f.R_s/f.d*(1+f.delta*f.rd -f.tau*f.delta*f.rtd);
end dsdT;
