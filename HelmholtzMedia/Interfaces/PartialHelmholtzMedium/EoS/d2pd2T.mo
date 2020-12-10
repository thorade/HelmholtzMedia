within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2pd2T "returns pressure derivative (d2p/dd2)@T=const"
  input EoS.HelmholtzDerivs f;
  output Der2PressureByDensity2 d2pd2T;

algorithm
  d2pd2T := f.T*f.R_s/f.d*(2*f.delta*f.rd + 4*f.delta*f.delta*f.rdd + f.delta*f.delta*f.delta*f.rddd);
end d2pd2T;
