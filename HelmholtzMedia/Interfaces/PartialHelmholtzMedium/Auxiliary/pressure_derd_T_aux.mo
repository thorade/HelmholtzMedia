within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function pressure_derd_T_aux "returns pressure derivative (dp/dd)@T=const (computed from the shared auxiliary properties, see pressure_derd_T)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerPressureByDensity dpdT;
algorithm
  if state.phase==1 then
    dpdT := EoS.dpdT(aux.f);
  elseif state.phase==2 then
    dpdT := Modelica.Constants.small; // zero
  else
    assert(false, "pressure_derd_T error: phase is neither 1 nor 2");
  end if;
end pressure_derd_T_aux;
