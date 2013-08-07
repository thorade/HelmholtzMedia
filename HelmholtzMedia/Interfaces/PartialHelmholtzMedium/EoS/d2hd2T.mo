within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function d2hd2T "returns enthalpy derivative (d2h/dd2)@T=const"
  input EoS.HelmholtzDerivs f;
  output Der2EnthalpyByDensity2 d2hd2T;

algorithm
  d2hd2T := f.R*f.T/f.d^2*((f.tau*f.delta*f.delta*f.rtdd) + (2*f.delta*f.delta*f.rdd + f.delta*f.delta*f.delta*f.rddd));
end d2hd2T;
