within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2ud2T "returns internal energy derivative (d2u/dd2)@T=const"
  input EoS.HelmholtzDerivs f;
  output Der2EnergyByDensity2 d2ud2T;

algorithm
  d2ud2T := f.R*f.T/f.d^2*(f.tau*f.delta*f.delta*f.rtdd);
end d2ud2T;
