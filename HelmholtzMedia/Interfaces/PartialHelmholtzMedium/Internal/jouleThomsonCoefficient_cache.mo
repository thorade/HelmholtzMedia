within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function jouleThomsonCoefficient_cache "returns Joule-Thomson-Coefficient (dT/dp)@h=const (computed from the shared auxiliary properties, see jouleThomsonCoefficient)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output DerTemperatureByPressure mu;
algorithm
  if state.phase==1 then
    mu := 1.0/(EoS.dpTd(cache.f)-EoS.dpdT(cache.f)*EoS.dhTd(cache.f)/EoS.dhdT(cache.f));
  elseif state.phase==2 then
    mu := saturationTemperature_derp_sat(sat=setSat_T(T=state.T));
  else
    assert(false, "jouleThomsonCoefficient error: phase is neither 1 nor 2");
  end if;
end jouleThomsonCoefficient_cache;
