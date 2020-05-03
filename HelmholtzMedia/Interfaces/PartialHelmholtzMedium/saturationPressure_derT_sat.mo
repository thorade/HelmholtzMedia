within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function saturationPressure_derT_sat "returns (dp/dT)@sat"

input SaturationProperties sat;
output DerPressureByTemperature dpT;

algorithm
  // vap entropy is larger than liq entropy, numerator is positive
  // vap volume is larger than liq volume, denominator is positive
  // Clausius-Clapeyron equation
  dpT := (sat.vap.s-sat.liq.s)/(1.0/sat.vap.d-1.0/sat.liq.d);
annotation(Inline = true);
end saturationPressure_derT_sat;
