within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derd_T "returns pressure derivative (dp/dd)@T=const"
  input ThermodynamicState state;
  output DerPressureByDensity dpdT;
algorithm
  dpdT := Internal.pressure_derd_T_cache(state=state, cache=Internal.setHelmholtzCache(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end pressure_derd_T;
