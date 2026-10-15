within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function dDewDensity_dPressure_cache "Return dew point density derivative (worker of dDewDensity_dPressure, with the derivatives at the saturated vapour density)"
  extends Modelica.Icons.Function;
  input Temperature T "saturation temperature, sat.Tsat (for the derivative annotation)";
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at sat.vap";
  output DerDensityByPressure ddvdp;
protected
  DerDensityByPressure ddpT = 1.0/EoS.dpdT(f);
  DerDensityByTemperature ddTp = -EoS.dpTd(f)/EoS.dpdT(f);
  DerTemperatureByPressure dTp = saturationTemperature_derp_sat(sat=sat);
algorithm
  ddvdp := ddpT + ddTp*dTp;
  annotation (derivative(noDerivative=sat, noDerivative=f)=dDewDensity_dPressure_cache_der);
end dDewDensity_dPressure_cache;
