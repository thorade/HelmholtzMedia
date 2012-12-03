within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2uT2d "returns internal energy derivative (d2u/dT2)@d=const"
  input EoS.HelmholtzDerivs f;
  output Der2EnergyByTemperature2 d2uT2d;

algorithm
  d2uT2d := f.R/f.T*(f.tau*f.tau*f.tau*(f.ittt + f.rttt) +2*f.tau*f.tau*(f.itt + f.rtt));
end d2uT2d;
