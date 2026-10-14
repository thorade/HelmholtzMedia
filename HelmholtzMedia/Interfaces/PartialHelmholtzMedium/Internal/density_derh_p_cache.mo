within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derh_p_cache "returns density derivative (dd/dh)@p=const (computed from the shared auxiliary properties, see density_derh_p)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output DerDensityByEnthalpy ddhp;
algorithm
  if state.phase==1 then
    ddhp := 1.0/(EoS.dhdT(cache.f) - EoS.dhTd(cache.f)*EoS.dpdT(cache.f)/EoS.dpTd(cache.f));
  elseif state.phase==2 then
    // dvhp = (v"-v')/(h"-h')
    // ddhp = -d^2 * dvhp
    ddhp := -state.d^2*(1/cache.sat.liq.d-1/cache.sat.vap.d)/(cache.sat.liq.h-cache.sat.vap.h);
  else
    assert(false, "density_derh_p error: phase is neither 1 nor 2");
  end if;
end density_derh_p_cache;
