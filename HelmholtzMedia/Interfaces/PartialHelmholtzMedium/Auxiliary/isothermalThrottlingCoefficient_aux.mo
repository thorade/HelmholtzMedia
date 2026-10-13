within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function isothermalThrottlingCoefficient_aux "returns isothermal throttling coefficient (dh/dp)@T=const (computed from the shared auxiliary properties, see isothermalThrottlingCoefficient)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerEnthalpyByPressure delta_T;
algorithm
  if state.phase==1 then
    delta_T := EoS.dhdT(aux.f)/EoS.dpdT(aux.f);
  elseif state.phase==2 then
    assert(false, "isothermalThrottlingCoefficient warning: in the two-phase region delta_T is infinite", level=AssertionLevel.warning);
    delta_T := Modelica.Constants.inf;
  else
    assert(false, "isothermalThrottlingCoefficient error: phase is neither 1 nor 2");
  end if;
end isothermalThrottlingCoefficient_aux;
