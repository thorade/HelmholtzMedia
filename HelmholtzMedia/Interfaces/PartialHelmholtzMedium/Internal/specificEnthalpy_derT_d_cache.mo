within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_derT_d_cache "returns enthalpy derivative (dh/dT)@d=const (computed from the shared auxiliary properties, see specificEnthalpy_derT_d)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerEnthalpyByTemperature dhTd;
protected
  MassFraction x "vapour quality";
  DerPressureByTemperature dpT;
  DerEnthalpyByTemperature dhT_liq;
  DerEnthalpyByTemperature dhT_vap;
  DerDensityByTemperature ddT_liq;
  DerDensityByTemperature ddT_vap;
  DerVolumeByTemperature dvT_liq;
  DerVolumeByTemperature dvT_vap;
  DerFractionByTemperature dxTv;
algorithm
  if phase==1 then
    dhTd := EoS.dhTd(cache.f);
  elseif phase==2 then
    x := (1/d - 1/cache.sat.liq.d)/(1/cache.sat.vap.d - 1/cache.sat.liq.d);
    dpT := (cache.sat.vap.s-cache.sat.liq.s)/(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d);
    dhT_liq := EoS.dhTd(cache.fl)-EoS.dhdT(cache.fl)*EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl) + EoS.dhdT(cache.fl)/EoS.dpdT(cache.fl)*dpT;
    dhT_vap := EoS.dhTd(cache.fv)-EoS.dhdT(cache.fv)*EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv) + EoS.dhdT(cache.fv)/EoS.dpdT(cache.fv)*dpT;
    ddT_liq := -EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl) + 1.0/EoS.dpdT(cache.fl)*dpT;
    ddT_vap := -EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv) + 1.0/EoS.dpdT(cache.fv)*dpT;
    dvT_liq := -1/cache.sat.liq.d^2 * ddT_liq;
    dvT_vap := -1/cache.sat.vap.d^2 * ddT_vap;
    dxTv :=(x*dvT_vap + (1 - x)*dvT_liq)/(1/cache.sat.liq.d - 1/cache.sat.vap.d);
    dhTd := dhT_liq + dxTv*(cache.sat.vap.h-cache.sat.liq.h) + x*(dhT_vap-dhT_liq);
  else
    assert(false, "specificEnthalpy_derT_d error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=specificEnthalpy_derT_d_cache_der);
end specificEnthalpy_derT_d_cache;
