within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derp_T_cache "returns density derivative (dd/dp)@T=const (computed from the shared auxiliary properties, see density_derp_T)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output DerDensityByPressure ddpT;
algorithm
  if state.phase==1 then
    ddpT := 1.0/EoS.dpdT(cache.f);
  elseif state.phase==2 then
    ddpT := Modelica.Constants.inf; // divide by zero
  end if;
end density_derp_T_cache;
