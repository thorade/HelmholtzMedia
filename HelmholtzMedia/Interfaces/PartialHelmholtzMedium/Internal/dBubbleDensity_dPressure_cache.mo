within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function dBubbleDensity_dPressure_cache "Return bubble point density derivative (worker of dBubbleDensity_dPressure, with the derivatives at the saturated liquid density)"
  extends Modelica.Icons.Function;
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at sat.liq";
  output DerDensityByPressure ddldp;
protected
  DerDensityByPressure ddpT = 1.0/EoS.dpdT(f);
  DerDensityByTemperature ddTp = -EoS.dpTd(f)/EoS.dpdT(f);
  DerTemperatureByPressure dTp = saturationTemperature_derp_sat(sat=sat);
algorithm
  ddldp := ddpT + ddTp*dTp;
end dBubbleDensity_dPressure_cache;
