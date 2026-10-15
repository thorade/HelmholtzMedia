within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function dBubbleEnthalpy_dPressure_cache "Return bubble point enthalpy derivative (worker of dBubbleEnthalpy_dPressure, with the derivatives at the saturated liquid density)"
  extends Modelica.Icons.Function;
  input Temperature T "saturation temperature, sat.Tsat (for the derivative annotation)";
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at sat.liq";
  output DerEnthalpyByPressure dhldp;
protected
  DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
  DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
  DerTemperatureByPressure dTp = saturationTemperature_derp_sat(sat=sat);
algorithm
  dhldp := dhpT + dhTp*dTp;
  annotation (derivative(noDerivative=sat, noDerivative=f)=dBubbleEnthalpy_dPressure_cache_der);
end dBubbleEnthalpy_dPressure_cache;
