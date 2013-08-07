within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dsTd "returns entropy derivative (ds/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEntropyByTemperature dsTd;

algorithm
  dsTd := f.R/f.T*(-f.tau*f.tau*(f.itt+f.rtt));
end dsTd;
