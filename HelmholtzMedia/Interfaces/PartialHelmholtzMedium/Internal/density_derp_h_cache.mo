within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derp_h_cache "returns density derivative (dd/dp)@h=const (computed from the shared auxiliary properties, see density_derp_h)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
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
    ddph := 1.0/(EoS.dpdT(cache.f) - EoS.dpTd(cache.f)*EoS.dhdT(cache.f)/EoS.dhTd(cache.f));
  elseif state.phase==2 then
    x := vapourQuality_sat(state=state, sat=cache.sat);
    dTp := saturationTemperature_derp_sat(sat=cache.sat);
    dhp_liq := EoS.dhdT(cache.fl)/EoS.dpdT(cache.fl) + (EoS.dhTd(cache.fl)-EoS.dhdT(cache.fl)*EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl))*dTp;
    dhp_vap := EoS.dhdT(cache.fv)/EoS.dpdT(cache.fv) + (EoS.dhTd(cache.fv)-EoS.dhdT(cache.fv)*EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv))*dTp;
    ddp_liq := 1.0/EoS.dpdT(cache.fl) - EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl)*dTp;
    ddp_vap := 1.0/EoS.dpdT(cache.fv) - EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv)*dTp;
    dvp_liq := -1.0/cache.sat.liq.d^2 * ddp_liq;
    dvp_vap := -1.0/cache.sat.vap.d^2 * ddp_vap;
    dxph :=(x*dhp_vap + (1 - x)*dhp_liq)/(cache.sat.liq.h - cache.sat.vap.h);
    dvph := dvp_liq + dxph*(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d) + x*(dvp_vap-dvp_liq);
    ddph := -state.d*state.d*dvph;
  end if;
end density_derp_h_cache;
