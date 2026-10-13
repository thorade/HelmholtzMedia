within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derT_d "returns pressure derivative (dp/dT)@d=const"
  input ThermodynamicState state;
  output DerPressureByTemperature dpTd;
algorithm
  dpTd := Internal.pressure_derT_d_cache(state=state, cache=Internal.setHelmholtzCache(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end pressure_derT_d;
