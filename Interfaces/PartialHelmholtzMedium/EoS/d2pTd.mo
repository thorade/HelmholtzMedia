within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2pTd "returns pressure derivative (d2p/dTdd)"
  input EoS.HelmholtzDerivs f;
  output Der2PressureByTemperatureDensity d2pTd;

algorithm
  d2pTd := f.R*(1 + 2*f.delta*f.rd + f.delta*f.delta*f.rdd - 2*f.delta*f.tau*f.rtd - f.tau*f.delta*f.delta*f.rtdd);
end d2pTd;
