within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2hT2d "returns enthalpy derivative (d2h/dT2)@d=const"
  input EoS.HelmholtzDerivs f;
  output Der2EnthalpyByTemperature2 d2hT2d;

algorithm
  d2hT2d := f.R/f.T*(f.tau*f.tau*f.tau*(f.ittt + f.rttt) + 2*f.tau*f.tau*(f.itt + f.rtt) + f.tau^2*f.delta*f.rttd);
end d2hT2d;
