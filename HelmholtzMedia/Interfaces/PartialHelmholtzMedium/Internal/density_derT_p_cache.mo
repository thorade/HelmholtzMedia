within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derT_p_cache "returns density derivative (dd/dT)@p=const (computed from the shared auxiliary properties, see density_derT_p)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output DerDensityByTemperature ddTp;
algorithm
  if state.phase==1 then
    ddTp := -EoS.dpTd(cache.f)/EoS.dpdT(cache.f);
  elseif state.phase==2 then
    ddTp := Modelica.Constants.inf; // divide by zero
  else
    assert(false, "density_derT_p error: phase is neither 1 nor 2");
  end if;
end density_derT_p_cache;
