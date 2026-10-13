within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function pressure_derT_d_cache "returns pressure derivative (dp/dT)@d=const (computed from the shared auxiliary properties, see pressure_derT_d)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output DerPressureByTemperature dpTd;
algorithm
  if state.phase==1 then
    dpTd := EoS.dpTd(cache.f);
  elseif state.phase==2 then
    dpTd := (cache.sat.vap.s-cache.sat.liq.s)/(1.0/cache.sat.vap.d-1.0/cache.sat.liq.d);
  else
    assert(false, "pressure_derT_d error: phase is neither 1 nor 2");
  end if;
end pressure_derT_d_cache;
