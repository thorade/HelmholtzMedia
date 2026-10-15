within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derT_p_cache "returns density derivative (dd/dT)@p=const (computed from the shared auxiliary properties, see density_derT_p)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerDensityByTemperature ddTp;
algorithm
  if phase==1 then
    ddTp := -EoS.dpTd(cache.f)/EoS.dpdT(cache.f);
  elseif phase==2 then
    ddTp := Modelica.Constants.inf; // divide by zero
  else
    assert(false, "density_derT_p error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=density_derT_p_cache_der);
end density_derT_p_cache;
