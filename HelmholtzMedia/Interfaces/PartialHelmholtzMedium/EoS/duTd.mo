within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function duTd "returns internal energy derivative (du/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByTemperature duTd;

algorithm
  duTd := f.R*(-f.tau*f.tau*(f.itt + f.rtt));
end duTd;
