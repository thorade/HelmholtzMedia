within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificHeatCapacityCp_cache "returns the isobaric specific heat capcacity (computed from the shared auxiliary properties, see specificHeatCapacityCp)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output SpecificHeatCapacity cp;
algorithm
  if phase==1 then
    cp := EoS.dhTd(cache.f) - EoS.dhdT(cache.f)*EoS.dpTd(cache.f)/EoS.dpdT(cache.f);
  elseif phase==2 then
    assert(false, "specificHeatCapacityCp warning: in the two-phase region cp is infinite", level=AssertionLevel.warning);
    cp := Modelica.Constants.inf;
  else
    assert(false, "specificHeatCapacityCp error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=specificHeatCapacityCp_cache_der);
end specificHeatCapacityCp_cache;
