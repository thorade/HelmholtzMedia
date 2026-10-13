within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function density_derT_p_aux "returns density derivative (dd/dT)@p=const (computed from the shared auxiliary properties, see density_derT_p)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerDensityByTemperature ddTp;
algorithm
  if state.phase==1 then
    ddTp := -EoS.dpTd(aux.f)/EoS.dpdT(aux.f);
  elseif state.phase==2 then
    ddTp := Modelica.Constants.inf; // divide by zero
  end if;
end density_derT_p_aux;
