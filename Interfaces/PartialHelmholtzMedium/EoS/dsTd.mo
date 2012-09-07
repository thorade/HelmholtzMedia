within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dsTd "returns pressure derivative (ds/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEntropyByTemperature dsTd;

algorithm
  dsTd := f.R/f.T*(-f.tau^2*(f.itt+f.rtt));
end dsTd;
