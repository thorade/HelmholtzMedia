within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function density_derp_T_aux "returns density derivative (dd/dp)@T=const (computed from the shared auxiliary properties, see density_derp_T)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerDensityByPressure ddpT;
algorithm
  if state.phase==1 then
    ddpT := 1.0/EoS.dpdT(aux.f);
  elseif state.phase==2 then
    ddpT := Modelica.Constants.inf; // divide by zero
  end if;
end density_derp_T_aux;
