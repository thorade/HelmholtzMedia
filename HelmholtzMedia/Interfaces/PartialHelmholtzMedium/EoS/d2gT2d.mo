within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2gT2d "returns Gibbs energy derivative (d2g/dT2)@d=const"
  input EoS.HelmholtzDerivs f;
  output Der2EnergyByTemperature2 d2gT2d;

algorithm
  d2gT2d := f.R/f.T*(f.tau*f.tau*(f.itt+f.rtt) + (f.tau^2*f.delta*f.rttd));
end d2gT2d;
