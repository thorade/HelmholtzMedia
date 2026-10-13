within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isothermalThrottlingCoefficient "returns isothermal throttling coefficient (dh/dp)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByPressure delta_T;
algorithm
  delta_T := isothermalThrottlingCoefficient_aux(state=state, aux=setAuxiliaryProperties(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end isothermalThrottlingCoefficient;
