within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function isothermalCompressibility_cache "returns -1/v*(dv/dp)@T=const (computed from the shared auxiliary properties, see isothermalCompressibility)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output Modelica.Units.SI.IsothermalCompressibility kappa;
algorithm
  if phase==1 then
    kappa := 1.0/d/EoS.dpdT(cache.f);
  elseif phase==2 then
    assert(false, "isothermalCompressibility warning: in the two-phase region kappa is infinite", level=AssertionLevel.warning);
    kappa := Modelica.Constants.inf;
  else
    assert(false, "isothermalCompressibility error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=isothermalCompressibility_cache_der);
end isothermalCompressibility_cache;
