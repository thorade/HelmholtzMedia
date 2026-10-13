within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function density_derh_p_aux "returns density derivative (dd/dh)@p=const (computed from the shared auxiliary properties, see density_derh_p)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerDensityByEnthalpy ddhp;
algorithm
  if state.phase==1 then
    ddhp := 1.0/(EoS.dhdT(aux.f) - EoS.dhTd(aux.f)*EoS.dpdT(aux.f)/EoS.dpTd(aux.f));
  elseif state.phase==2 then
    // dvhp = (v"-v')/(h"-h')
    // ddhp = -d^2 * dvhp
    ddhp := -state.d^2*(1/aux.sat.liq.d-1/aux.sat.vap.d)/(aux.sat.liq.h-aux.sat.vap.h);
  end if;
end density_derh_p_aux;
