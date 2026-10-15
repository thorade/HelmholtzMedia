within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function pressure_derd_T_cache "returns pressure derivative (dp/dd)@T=const (computed from the shared auxiliary properties, see pressure_derd_T)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerPressureByDensity dpdT;
algorithm
  if phase==1 then
    dpdT := EoS.dpdT(cache.f);
  elseif phase==2 then
    dpdT := Modelica.Constants.small; // zero
  else
    assert(false, "pressure_derd_T error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=pressure_derd_T_cache_der);
end pressure_derd_T_cache;
