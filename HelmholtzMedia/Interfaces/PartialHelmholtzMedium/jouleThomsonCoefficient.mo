within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function jouleThomsonCoefficient "returns Joule-Thomson-Coefficient (dT/dp)@h=const"
  input ThermodynamicState state;
  output DerTemperatureByPressure mu;
algorithm
  mu := jouleThomsonCoefficient_aux(state=state, aux=setAuxiliaryProperties(d=state.d, T=state.T, phase=state.phase));
annotation (Inline=true);
end jouleThomsonCoefficient;
