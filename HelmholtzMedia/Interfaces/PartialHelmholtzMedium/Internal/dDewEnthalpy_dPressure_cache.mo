within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function dDewEnthalpy_dPressure_cache "Return dew point enthalpy derivative (worker of dDewEnthalpy_dPressure, with the derivatives at the saturated vapour density)"
  extends Modelica.Icons.Function;
  input Temperature T "saturation temperature, sat.Tsat (for the derivative annotation)";
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at sat.vap";
  output DerEnthalpyByPressure dhvdp;
protected
  DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
  DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
  DerTemperatureByPressure dTp = saturationTemperature_derp_sat(sat=sat);
algorithm
  dhvdp := dhpT + dhTp*dTp;
  annotation (derivative(noDerivative=sat, noDerivative=f)=dDewEnthalpy_dPressure_cache_der);
end dDewEnthalpy_dPressure_cache;
