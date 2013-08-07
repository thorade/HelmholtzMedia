within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2sd2T "returns entropy derivative (d2s/dd2)@T=const"
  input EoS.HelmholtzDerivs f;
  output Der2EntropyByDensity2 d2sd2T;

algorithm
  d2sd2T := f.R/f.d^2*(1 - f.delta*f.delta*f.rdd + f.tau*f.delta*f.delta*f.rtdd);
end d2sd2T;
