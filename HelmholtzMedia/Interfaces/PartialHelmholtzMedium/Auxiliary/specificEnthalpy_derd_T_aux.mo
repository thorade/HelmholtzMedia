within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function specificEnthalpy_derd_T_aux "returns enthalpy derivative (dh/dd)@T=const (computed from the shared auxiliary properties, see specificEnthalpy_derd_T)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerEnthalpyByDensity dhdT;
algorithm
  if state.phase==1 then
    dhdT := EoS.dhdT(aux.f);
  elseif state.phase==2 then
    // dhvT = (h"-h')/(v"-v')
    // dhdT = -1/d^2 * dhvT
    dhdT := -1/state.d^2*(aux.sat.liq.h-aux.sat.vap.h)/(1/aux.sat.liq.d-1/aux.sat.vap.d);
  else
    assert(false, "specificEnthalpy_derd_T error: phase is neither 1 nor 2");
  end if;
end specificEnthalpy_derd_T_aux;
