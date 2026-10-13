within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function specificHeatCapacityCp_aux "returns the isobaric specific heat capcacity (computed from the shared auxiliary properties, see specificHeatCapacityCp)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output SpecificHeatCapacity cp;
algorithm
  if state.phase==1 then
    cp := EoS.dhTd(aux.f) - EoS.dhdT(aux.f)*EoS.dpTd(aux.f)/EoS.dpdT(aux.f);
  elseif state.phase==2 then
    assert(false, "specificHeatCapacityCp warning: in the two-phase region cp is infinite", level=AssertionLevel.warning);
    cp := Modelica.Constants.inf;
  end if;
end specificHeatCapacityCp_aux;
