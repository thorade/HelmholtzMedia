within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function velocityOfSound_cache "returns the speed or velocity of sound (computed from the shared auxiliary properties, see velocityOfSound)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
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
  if phase==1 then
    a := sqrt(EoS.dpdT(cache.f)-EoS.dpTd(cache.f)*EoS.dsdT(cache.f)/EoS.dsTd(cache.f));
  elseif phase==2 then
    x := if d <= cache.sat.vap.d then 1.0 elseif d >= cache.sat.liq.d then 0.0
         else (1.0/d - 1.0/cache.sat.liq.d)/(1.0/cache.sat.vap.d - 1.0/cache.sat.liq.d);
    dTp := saturationTemperature_derp_sat(sat=cache.sat);
    dsp_liq := EoS.dsdT(cache.fl)/EoS.dpdT(cache.fl) + (EoS.dsTd(cache.fl)-EoS.dsdT(cache.fl)*EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl))*dTp;
    dsp_vap := EoS.dsdT(cache.fv)/EoS.dpdT(cache.fv) + (EoS.dsTd(cache.fv)-EoS.dsdT(cache.fv)*EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv))*dTp;
    ddp_liq := 1.0/EoS.dpdT(cache.fl) - EoS.dpTd(cache.fl)/EoS.dpdT(cache.fl)*dTp;
    ddp_vap := 1.0/EoS.dpdT(cache.fv) - EoS.dpTd(cache.fv)/EoS.dpdT(cache.fv)*dTp;
    dvp_liq := -1.0/cache.sat.liq.d^2 * ddp_liq;
    dvp_vap := -1.0/cache.sat.vap.d^2 * ddp_vap;
    dxps :=(x*dsp_vap + (1 - x)*dsp_liq)/(cache.sat.liq.s - cache.sat.vap.s);
    dvps := dvp_liq + dxps*(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d) + x*(dvp_vap-dvp_liq);
    a := sqrt(-1.0/(d*d*dvps));
  else
    assert(false, "velocityOfSound error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=velocityOfSound_cache_der);
end velocityOfSound_cache;
