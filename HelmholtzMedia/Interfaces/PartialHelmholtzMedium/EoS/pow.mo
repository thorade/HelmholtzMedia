within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function pow "power a^b"

  input Real a;
  input Real b;
  output Real y;

external "C" y = pow(a,b);

end pow;
