within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2pTv "returns pressure derivative (d2p/dT dv)"
  input EoS.HelmholtzDerivs f;
  output Der2PressureByTemperatureVolume d2pTv;

algorithm
  d2pTv := -f.d*f.d*f.R*(1 + 2*f.delta*f.rd + f.delta*f.delta*f.rdd - 2*f.delta*f.tau*f.rtd - f.tau*f.delta*f.delta*f.rtdd);
end d2pTv;
