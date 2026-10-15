within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derp_T_cache "returns density derivative (dd/dp)@T=const (computed from the shared auxiliary properties, see density_derp_T)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output DerDensityByPressure ddpT;
algorithm
  if phase==1 then
    ddpT := 1.0/EoS.dpdT(cache.f);
  elseif phase==2 then
    ddpT := Modelica.Constants.inf; // divide by zero
  else
    assert(false, "density_derp_T error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=density_derp_T_cache_der);
end density_derp_T_cache;
