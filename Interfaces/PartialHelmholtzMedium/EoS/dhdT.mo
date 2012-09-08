within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dhdT "returns enthalpy derivative (dh/dd)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerEnthalpyByDensity dhdT;

algorithm
  dhdT := f.T*f.R/f.d*(0+f.tau*f.delta*f.rtd + f.delta*f.rd + f.delta^2*f.rdd);
end dhdT;
