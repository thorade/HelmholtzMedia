within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derd_T "returns pressure derivative (dp/dd)@T=const"
  input ThermodynamicState state;
  output DerPressureByDensity dpdT;
algorithm
  dpdT := Auxiliary.pressure_derd_T_aux(state=state, aux=Auxiliary.setAuxiliaryProperties(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end pressure_derd_T;
