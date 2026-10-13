within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function velocityOfSound_cache "returns the speed or velocity of sound (computed from the shared auxiliary properties, see velocityOfSound)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output VelocityOfSound a;
protected
  MassFraction x "vapour quality";
  DerTemperatureByPressure dTp;
  DerEntropyByPressure dsp_liq;
  DerEntropyByPressure dsp_vap;
  DerDensityByPressure ddp_liq;
  DerDensityByPressure ddp_vap;
  DerVolumeByPressure dvp_liq;
  DerVolumeByPressure dvp_vap;
  DerVolumeByPressure dvps;
  DerFractionByPressure dxps;
algorithm
  if state.phase==1 then
    a := sqrt(EoS.dpdT(cache.f)-EoS.dpTd(cache.f)*EoS.dsdT(cache.f)/EoS.dsTd(cache.f));
  elseif state.phase==2 then
    x := vapourQuality_sat(state=state, sat=cache.sat);
    dTp := saturationTemperature_derp_sat(sat=cache.sat);
    dsp_liq := EoS.dsdT(cache.fl)/EoS.dpdT(cache.fl) + (EoS.dsTd(cache.fl)-EoS.dsdT(cache.fl)*EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl))*dTp;
    dsp_vap := EoS.dsdT(cache.fv)/EoS.dpdT(cache.fv) + (EoS.dsTd(cache.fv)-EoS.dsdT(cache.fv)*EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv))*dTp;
    ddp_liq := 1.0/EoS.dpdT(cache.fl) - EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl)*dTp;
    ddp_vap := 1.0/EoS.dpdT(cache.fv) - EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv)*dTp;
    dvp_liq := -1.0/cache.sat.liq.d^2 * ddp_liq;
    dvp_vap := -1.0/cache.sat.vap.d^2 * ddp_vap;
    dxps :=(x*dsp_vap + (1 - x)*dsp_liq)/(cache.sat.liq.s - cache.sat.vap.s);
    dvps := dvp_liq + dxps*(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d) + x*(dvp_vap-dvp_liq);
    a := sqrt(-1.0/(state.d*state.d*dvps));
  end if;
end velocityOfSound_cache;
