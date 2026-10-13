within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function isothermalCompressibility_aux "returns -1/v*(dv/dp)@T=const (computed from the shared auxiliary properties, see isothermalCompressibility)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output Modelica.Units.SI.IsothermalCompressibility kappa;
algorithm
  if state.phase==1 then
    kappa := 1.0/state.d/EoS.dpdT(aux.f);
  elseif state.phase==2 then
    assert(false, "isothermalCompressibility warning: in the two-phase region kappa is infinite", level=AssertionLevel.warning);
    kappa := Modelica.Constants.inf;
  end if;
end isothermalCompressibility_aux;
