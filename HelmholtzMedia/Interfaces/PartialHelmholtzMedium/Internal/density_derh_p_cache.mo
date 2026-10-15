within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derh_p_cache "returns density derivative (dd/dh)@p=const (computed from the shared auxiliary properties, see density_derh_p)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerDensityByEnthalpy ddhp;
algorithm
  if phase==1 then
    ddhp := 1.0/(EoS.dhdT(cache.f) - EoS.dhTd(cache.f)*EoS.dpdT(cache.f)/EoS.dpTd(cache.f));
  elseif phase==2 then
    // dvhp = (v"-v')/(h"-h')
    // ddhp = -d^2 * dvhp
    ddhp := -d^2*(1/cache.sat.liq.d-1/cache.sat.vap.d)/(cache.sat.liq.h-cache.sat.vap.h);
  else
    assert(false, "density_derh_p error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=density_derh_p_cache_der);
end density_derh_p_cache;
