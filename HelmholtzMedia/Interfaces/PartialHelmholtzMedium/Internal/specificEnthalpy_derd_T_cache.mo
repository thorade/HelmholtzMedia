within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_derd_T_cache "returns enthalpy derivative (dh/dd)@T=const (computed from the shared auxiliary properties, see specificEnthalpy_derd_T)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerEnthalpyByDensity dhdT;
algorithm
  if phase==1 then
    dhdT := EoS.dhdT(cache.f);
  elseif phase==2 then
    // dhvT = (h"-h')/(v"-v')
    // dhdT = -1/d^2 * dhvT
    dhdT := -1/d^2*(cache.sat.liq.h-cache.sat.vap.h)/(1/cache.sat.liq.d-1/cache.sat.vap.d);
  else
    assert(false, "specificEnthalpy_derd_T error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=specificEnthalpy_derd_T_cache_der);
end specificEnthalpy_derd_T_cache;
