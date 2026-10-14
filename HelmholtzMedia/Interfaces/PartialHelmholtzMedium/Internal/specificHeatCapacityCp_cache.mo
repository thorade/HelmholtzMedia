within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificHeatCapacityCp_cache "returns the isobaric specific heat capcacity (computed from the shared auxiliary properties, see specificHeatCapacityCp)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output SpecificHeatCapacity cp;
algorithm
  if state.phase==1 then
    cp := EoS.dhTd(cache.f) - EoS.dhdT(cache.f)*EoS.dpTd(cache.f)/EoS.dpdT(cache.f);
  elseif state.phase==2 then
    assert(false, "specificHeatCapacityCp warning: in the two-phase region cp is infinite", level=AssertionLevel.warning);
    cp := Modelica.Constants.inf;
  else
    assert(false, "specificHeatCapacityCp error: phase is neither 1 nor 2");
  end if;
end specificHeatCapacityCp_cache;
