within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2pT2d "returns pressure derivative (d2p/dT2)@d=const"
  input EoS.HelmholtzDerivs f;
  output Der2PressureByTemperature2 d2pT2d;

algorithm
  d2pT2d := f.d*f.R/f.T*(f.tau^2*f.delta*f.rttd);
end d2pT2d;
