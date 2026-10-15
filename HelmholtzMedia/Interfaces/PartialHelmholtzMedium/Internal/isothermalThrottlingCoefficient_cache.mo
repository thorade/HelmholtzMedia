within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function isothermalThrottlingCoefficient_cache "returns isothermal throttling coefficient (dh/dp)@T=const (computed from the shared auxiliary properties, see isothermalThrottlingCoefficient)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerEnthalpyByPressure delta_T;
algorithm
  if phase==1 then
    delta_T := EoS.dhdT(cache.f)/EoS.dpdT(cache.f);
  elseif phase==2 then
    assert(false, "isothermalThrottlingCoefficient warning: in the two-phase region delta_T is infinite", level=AssertionLevel.warning);
    delta_T := Modelica.Constants.inf;
  else
    assert(false, "isothermalThrottlingCoefficient error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=isothermalThrottlingCoefficient_cache_der);
end isothermalThrottlingCoefficient_cache;
