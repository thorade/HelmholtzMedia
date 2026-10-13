within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function jouleThomsonCoefficient_aux "returns Joule-Thomson-Coefficient (dT/dp)@h=const (computed from the shared auxiliary properties, see jouleThomsonCoefficient)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerTemperatureByPressure mu;
algorithm
  if state.phase==1 then
    mu := 1.0/(EoS.dpTd(aux.f)-EoS.dpdT(aux.f)*EoS.dhTd(aux.f)/EoS.dhdT(aux.f));
  elseif state.phase==2 then
    mu := saturationTemperature_derp_sat(sat=setSat_T(T=state.T));
  else
    assert(false, "jouleThomsonCoefficient error: phase is neither 1 nor 2");
  end if;
end jouleThomsonCoefficient_aux;
