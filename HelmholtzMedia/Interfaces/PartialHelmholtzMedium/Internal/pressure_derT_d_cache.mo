within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function pressure_derT_d_cache "returns pressure derivative (dp/dT)@d=const (computed from the shared auxiliary properties, see pressure_derT_d)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerPressureByTemperature dpTd;
algorithm
  if phase==1 then
    dpTd := EoS.dpTd(cache.f);
  elseif phase==2 then
    dpTd := (cache.sat.vap.s-cache.sat.liq.s)/(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d);
  else
    assert(false, "pressure_derT_d error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=pressure_derT_d_cache_der);
end pressure_derT_d_cache;
