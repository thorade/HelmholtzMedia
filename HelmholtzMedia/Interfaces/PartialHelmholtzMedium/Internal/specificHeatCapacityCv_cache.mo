within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificHeatCapacityCv_cache "returns the isochoric specific heat capcacity (computed from the shared auxiliary properties, see specificHeatCapacityCv)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output SpecificHeatCapacity cv;
protected
  MassFraction x "vapour quality";
  DerPressureByTemperature dpT;
  DerEnergyByTemperature duT_liq;
  DerEnergyByTemperature duT_vap;
  DerDensityByTemperature ddT_liq;
  DerDensityByTemperature ddT_vap;
  DerVolumeByTemperature dvT_liq;
  DerVolumeByTemperature dvT_vap;
  DerFractionByTemperature dxTv;
algorithm
  if state.phase==1 then
    cv := EoS.duTd(cache.f);
  elseif state.phase==2 then
    x := vapourQuality_sat(state=state, sat=cache.sat);
    dpT := saturationPressure_derT_sat(sat=cache.sat);
    duT_liq := EoS.duTd(cache.fl)-EoS.dudT(cache.fl)*EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl) + EoS.dudT(cache.fl)/EoS.dpdT(cache.fl)*dpT;
    duT_vap := EoS.duTd(cache.fv)-EoS.dudT(cache.fv)*EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv) + EoS.dudT(cache.fv)/EoS.dpdT(cache.fv)*dpT;
    ddT_liq := -EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl) + 1.0/EoS.dpdT(cache.fl)*dpT;
    ddT_vap := -EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv) + 1.0/EoS.dpdT(cache.fv)*dpT;
    dvT_liq := -1/cache.sat.liq.d^2 * ddT_liq;
    dvT_vap := -1/cache.sat.vap.d^2 * ddT_vap;
    dxTv :=(x*dvT_vap + (1 - x)*dvT_liq)/(1/cache.sat.liq.d - 1/cache.sat.vap.d);
    cv := duT_liq + dxTv*(cache.sat.vap.u-cache.sat.liq.u) + x*(duT_vap-duT_liq);
  else
    assert(false, "specificHeatCapacityCv error: phase is neither 1 nor 2");
  end if;
end specificHeatCapacityCv_cache;
