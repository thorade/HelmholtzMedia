within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function density_derp_h_aux "returns density derivative (dd/dp)@h=const (computed from the shared auxiliary properties, see density_derp_h)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerDensityByPressure ddph;
protected
  MassFraction x "vapour quality";
  DerTemperatureByPressure dTp;
  DerEnthalpyByPressure dhp_liq;
  DerEnthalpyByPressure dhp_vap;
  DerDensityByPressure ddp_liq;
  DerDensityByPressure ddp_vap;
  DerVolumeByPressure dvp_liq;
  DerVolumeByPressure dvp_vap;
  DerVolumeByPressure dvph;
  DerFractionByPressure dxph;
algorithm
  if state.phase==1 then
    ddph := 1.0/(EoS.dpdT(aux.f) - EoS.dpTd(aux.f)*EoS.dhdT(aux.f)/EoS.dhTd(aux.f));
  elseif state.phase==2 then
    x := vapourQuality_sat(state=state, sat=aux.sat);
    dTp := saturationTemperature_derp_sat(sat=aux.sat);
    dhp_liq := EoS.dhdT(aux.fl)/EoS.dpdT(aux.fl) + (EoS.dhTd(aux.fl)-EoS.dhdT(aux.fl)*EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl))*dTp;
    dhp_vap := EoS.dhdT(aux.fv)/EoS.dpdT(aux.fv) + (EoS.dhTd(aux.fv)-EoS.dhdT(aux.fv)*EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv))*dTp;
    ddp_liq := 1.0/EoS.dpdT(aux.fl) - EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl)*dTp;
    ddp_vap := 1.0/EoS.dpdT(aux.fv) - EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv)*dTp;
    dvp_liq := -1.0/aux.sat.liq.d^2 * ddp_liq;
    dvp_vap := -1.0/aux.sat.vap.d^2 * ddp_vap;
    dxph :=(x*dhp_vap + (1 - x)*dhp_liq)/(aux.sat.liq.h - aux.sat.vap.h);
    dvph := dvp_liq + dxph*(1.0/aux.sat.vap.d-1.0/aux.sat.liq.d) + x*(dvp_vap-dvp_liq);
    ddph := -state.d*state.d*dvph;
  end if;
end density_derp_h_aux;
