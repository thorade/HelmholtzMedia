within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2hTd "returns enthalpy derivative (d2h/dT dd)"
  input EoS.HelmholtzDerivs f;
  output Der2EnthalpyByTemperatureDensity d2hTd;

algorithm
  d2hTd := f.R_s/f.d*((-f.tau*f.tau*f.delta*f.rttd) + (f.delta*f.delta*f.rdd + f.delta*f.rd - f.tau*f.delta*f.delta*f.rtdd - f.tau*f.delta*f.rtd));
end d2hTd;
