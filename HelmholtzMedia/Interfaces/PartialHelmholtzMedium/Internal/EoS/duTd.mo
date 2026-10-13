within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.EoS;
function duTd "returns internal energy derivative (du/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByTemperature duTd;

algorithm
  duTd := f.R_s*(-f.tau*f.tau*(f.itt + f.rtt));
end duTd;
