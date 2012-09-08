within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function saturationPressure_derT "returns (dp/dT)@sat"

input Temperature T;
output DerPressureByTemperature dpT;

protected
  SaturationProperties sat=setSat_T(T=T);

algorithm
// Clausius-Clapeyron equation
   dpT := (sat.vap.s-sat.liq.s)/(1.0/sat.vap.d-1.0/sat.liq.d);
end saturationPressure_derT;
