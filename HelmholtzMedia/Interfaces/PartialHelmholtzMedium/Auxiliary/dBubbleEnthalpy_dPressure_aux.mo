within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function dBubbleEnthalpy_dPressure_aux "Return bubble point enthalpy derivative (worker of dBubbleEnthalpy_dPressure, with the derivatives at the saturated liquid density)"
  extends Modelica.Icons.Function;
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at sat.liq";
  output DerEnthalpyByPressure dhldp;
protected
  DerEnthalpyByPressure dhpT = EoS.dhdT(f)/EoS.dpdT(f);
  DerEnthalpyByTemperature dhTp = EoS.dhTd(f) - EoS.dhdT(f)*EoS.dpTd(f)/EoS.dpdT(f);
  DerTemperatureByPressure dTp = saturationTemperature_derp_sat(sat=sat);
algorithm
  dhldp := dhpT + dhTp*dTp;
end dBubbleEnthalpy_dPressure_aux;
