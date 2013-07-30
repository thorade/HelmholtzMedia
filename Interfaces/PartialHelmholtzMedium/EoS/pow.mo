within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function pow "power function with treatment of 0^0"
  input Real a "base";
  input Real b "exponent";
  output Real y;
algorithm
  y := if ((a<>0) or (b<>0)) then a^b else 1;
  annotation(Inline=true);
end pow;
