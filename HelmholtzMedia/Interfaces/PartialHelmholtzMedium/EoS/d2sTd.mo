within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2sTd "returns entropy derivative (d2s/dT dd)"
  input EoS.HelmholtzDerivs f;
  output Der2EntropyByTemperatureDensity d2sTd;

algorithm
  d2sTd := f.R_s/(f.T*f.d)*(-f.tau*f.tau*f.delta*f.rttd);
end d2sTd;
