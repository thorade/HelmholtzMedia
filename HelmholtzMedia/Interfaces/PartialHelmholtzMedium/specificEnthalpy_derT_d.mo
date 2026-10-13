within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derT_d "returns enthalpy derivative (dh/dT)@d=const"
  input ThermodynamicState state;
  output DerEnthalpyByTemperature dhTd;
algorithm
  dhTd := specificEnthalpy_derT_d_aux(state=state, aux=setAuxiliaryProperties(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end specificEnthalpy_derT_d;
