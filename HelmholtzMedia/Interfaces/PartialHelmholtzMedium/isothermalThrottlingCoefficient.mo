within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isothermalThrottlingCoefficient "returns isothermal throttling coefficient (dh/dp)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByPressure delta_T;
algorithm
  delta_T := Internal.isothermalThrottlingCoefficient_cache(d=state.d, T=state.T, phase=state.phase, cache=Internal.setHelmholtzCache(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end isothermalThrottlingCoefficient;
