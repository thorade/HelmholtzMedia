within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function saturationPressure_derT "returns (dp/dT)@sat"

input Temperature T;
output DerPressureByTemperature dpT;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;

algorithm
  dpT := if T< T_crit then saturationPressure_derT_sat(sat=setSat_T(T=T)) else pressure_derT_d(state=setState_dT(T=T, d=d_crit));
annotation(Inline = true);
end saturationPressure_derT;
