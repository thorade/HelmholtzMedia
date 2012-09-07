within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function duTd "returns pressure derivative (du/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByTemperature duTd;

algorithm
  duTd := f.R*(-f.tau^2*(f.itt + f.rtt));
end duTd;
