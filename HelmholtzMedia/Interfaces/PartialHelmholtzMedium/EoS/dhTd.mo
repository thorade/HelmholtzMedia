within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dhTd "returns enthalpy derivative (dh/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEnthalpyByTemperature dhTd;

algorithm
  dhTd := f.R_s*(1 - f.tau*f.tau*(f.itt+f.rtt) + f.delta*f.rd - f.tau*f.delta*f.rtd);
end dhTd;
