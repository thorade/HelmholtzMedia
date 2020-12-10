within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2pv2T "returns pressure derivative (d2p/dv2)@T=const"
  input EoS.HelmholtzDerivs f;
  output Der2PressureByVolume2 d2pv2T;

algorithm
  d2pv2T := f.d*f.d*f.d*f.T*f.R_s*(2 + 6*f.delta*f.rd + 6*f.delta*f.delta*f.rdd + f.delta*f.delta*f.delta*f.rddd);
end d2pv2T;
