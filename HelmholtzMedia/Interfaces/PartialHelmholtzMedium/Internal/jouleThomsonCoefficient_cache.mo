within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function jouleThomsonCoefficient_cache "returns Joule-Thomson-Coefficient (dT/dp)@h=const (computed from the shared auxiliary properties, see jouleThomsonCoefficient)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerTemperatureByPressure mu;
algorithm
  if phase==1 then
    mu := 1.0/(EoS.dpTd(cache.f)-EoS.dpdT(cache.f)*EoS.dhTd(cache.f)/EoS.dhdT(cache.f));
  elseif phase==2 then
    mu := saturationTemperature_derp_sat(sat=cache.sat);
  else
    assert(false, "jouleThomsonCoefficient error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=jouleThomsonCoefficient_cache_der);
end jouleThomsonCoefficient_cache;
