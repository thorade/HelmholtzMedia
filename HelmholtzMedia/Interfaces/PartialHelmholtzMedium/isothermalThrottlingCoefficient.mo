within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isothermalThrottlingCoefficient "returns isothermal throttling coefficient (dh/dp)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByPressure delta_T;
algorithm
  delta_T := Auxiliary.isothermalThrottlingCoefficient_aux(state=state, aux=Auxiliary.setAuxiliaryProperties(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end isothermalThrottlingCoefficient;
