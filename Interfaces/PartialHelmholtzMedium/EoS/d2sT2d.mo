within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2sT2d "returns entropy derivative (d2s/dT2)@d=const"
  input EoS.HelmholtzDerivs f;
  output Der2EntropyByTemperature2 d2sT2d;

algorithm
  d2sT2d := f.R/f.T^2*(f.tau*f.tau*f.tau*(f.ittt + f.rttt) + 2*f.tau*f.tau*(f.itt + f.rtt) + f.tau*f.tau*(f.itt+f.rtt));
end d2sT2d;
