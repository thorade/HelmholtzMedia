within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2gd2T "returns Gibbs energy derivative (d2g/dd2)@T=const"
  input EoS.HelmholtzDerivs f;
  output Der2EnergyByDensity2 d2gd2T;

algorithm
  d2gd2T := f.T*f.R/f.d^2*(-1 + (3*f.delta*f.delta*f.rdd + f.delta*f.delta*f.delta*f.rddd));
end d2gd2T;
