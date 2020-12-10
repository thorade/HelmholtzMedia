within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dudT "returns internal energy derivative (du/dd)@T=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByDensity dudT;

algorithm
  dudT := f.R_s*f.T/f.d*(f.tau*f.delta*f.rtd);
end dudT;
