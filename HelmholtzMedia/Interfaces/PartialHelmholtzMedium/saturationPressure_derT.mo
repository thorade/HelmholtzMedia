within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function saturationPressure_derT "returns (dp/dT)@sat"

input Temperature T;
output DerPressureByTemperature dpT;

algorithm
  dpT := saturationPressure_derT_sat(sat=setSat_T(T=T));
annotation(Inline = true);
end saturationPressure_derT;
