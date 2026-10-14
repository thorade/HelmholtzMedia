within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derp_h_cache "returns density derivative (dd/dp)@h=const (computed from the shared auxiliary properties, see density_derp_h)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
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
  if phase==1 then
    ddph := 1.0/(EoS.dpdT(cache.f) - EoS.dpTd(cache.f)*EoS.dhdT(cache.f)/EoS.dhTd(cache.f));
  elseif phase==2 then
    x := if d <= cache.sat.vap.d then 1.0 elseif d >= cache.sat.liq.d then 0.0
         else (1.0/d - 1.0/cache.sat.liq.d)/(1.0/cache.sat.vap.d - 1.0/cache.sat.liq.d);
    dTp := saturationTemperature_derp_sat(sat=cache.sat);
    dhp_liq := EoS.dhdT(cache.fl)/EoS.dpdT(cache.fl) + (EoS.dhTd(cache.fl)-EoS.dhdT(cache.fl)*EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl))*dTp;
    dhp_vap := EoS.dhdT(cache.fv)/EoS.dpdT(cache.fv) + (EoS.dhTd(cache.fv)-EoS.dhdT(cache.fv)*EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv))*dTp;
    ddp_liq := 1.0/EoS.dpdT(cache.fl) - EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl)*dTp;
    ddp_vap := 1.0/EoS.dpdT(cache.fv) - EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv)*dTp;
    dvp_liq := -1.0/cache.sat.liq.d^2 * ddp_liq;
    dvp_vap := -1.0/cache.sat.vap.d^2 * ddp_vap;
    dxph :=(x*dhp_vap + (1 - x)*dhp_liq)/(cache.sat.liq.h - cache.sat.vap.h);
    dvph := dvp_liq + dxph*(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d) + x*(dvp_vap-dvp_liq);
    ddph := -d*d*dvph;
  else
    assert(false, "density_derp_h error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=density_derp_h_cache_der);
end density_derp_h_cache;
