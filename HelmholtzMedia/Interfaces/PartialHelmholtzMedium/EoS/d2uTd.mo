within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2uTd "returns internal energy derivative (d2u/dT dd)"
  input EoS.HelmholtzDerivs f;
  output Der2EnergyByTemperatureDensity d2uTd;

algorithm
  d2uTd := f.R_s/f.d*(-f.tau*f.tau*f.delta*f.rttd);
end d2uTd;
