within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derT_d_aux "returns pressure derivative (dp/dT)@d=const (computed from the shared auxiliary properties, see pressure_derT_d)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerPressureByTemperature dpTd;
algorithm
  if state.phase==1 then
    dpTd := EoS.dpTd(aux.f);
  elseif state.phase==2 then
    dpTd := (aux.sat.vap.s-aux.sat.liq.s)/(1.0/aux.sat.vap.d-1.0/aux.sat.liq.d);
  else
    assert(false, "pressure_derT_d error: phase is neither 1 nor 2");
  end if;
end pressure_derT_d_aux;
