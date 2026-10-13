within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derd_T "returns enthalpy derivative (dh/dd)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByDensity dhdT;
algorithm
  dhdT := specificEnthalpy_derd_T_aux(state=state, aux=setAuxiliaryProperties(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end specificEnthalpy_derd_T;
