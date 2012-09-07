within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dhTd "returns pressure derivative (dh/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEnthalpyByTemperature dhTd;

algorithm
  dhTd := f.R*(1 - f.tau^2*(f.itt+f.rtt) + f.delta*f.rd - f.tau*f.delta*f.rtd);
end dhTd;
