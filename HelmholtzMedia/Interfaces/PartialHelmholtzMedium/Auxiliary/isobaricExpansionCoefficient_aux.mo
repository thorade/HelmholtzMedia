within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function isobaricExpansionCoefficient_aux "returns 1/v*(dv/dT)@p=const (computed from the shared auxiliary properties, see isobaricExpansionCoefficient)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output IsobaricExpansionCoefficient beta;
algorithm
  if state.phase==1 then
    beta := 1.0/state.d*EoS.dpTd(aux.f)/EoS.dpdT(aux.f);
  elseif state.phase==2 then
    assert(false, "isobaricExpansionCoefficient warning: in the two-phase region beta is zero", level=AssertionLevel.warning);
    beta := Modelica.Constants.small;
  end if;
end isobaricExpansionCoefficient_aux;
