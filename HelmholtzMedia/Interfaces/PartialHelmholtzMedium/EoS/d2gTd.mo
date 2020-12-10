within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2gTd "returns Gibbs energy derivative (d2g/dT dd)"
  input EoS.HelmholtzDerivs f;
  output Der2EnergyByTemperatureDensity d2gTd;

algorithm
  d2gTd := f.R_s/f.d*(1 + 2*f.delta*f.rd - 2*f.tau*f.delta*f.rtd + f.delta*f.delta*f.rdd - f.tau*f.delta*f.delta*f.rtdd);
end d2gTd;
