within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dpTd "returns pressure derivative (dp/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerPressureByTemperature dpTd;

algorithm
  dpTd := f.d*f.R*(1 + f.delta*f.rd - f.delta*f.tau*f.rtd);
end dpTd;
