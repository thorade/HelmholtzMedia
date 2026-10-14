within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function isobaricExpansionCoefficient_cache "returns 1/v*(dv/dT)@p=const (computed from the shared auxiliary properties, see isobaricExpansionCoefficient)"
  input ThermodynamicState state;
  input HelmholtzCache cache;
  output IsobaricExpansionCoefficient beta;
algorithm
  if state.phase==1 then
    beta := 1.0/state.d*EoS.dpTd(cache.f)/EoS.dpdT(cache.f);
  elseif state.phase==2 then
    assert(false, "isobaricExpansionCoefficient warning: in the two-phase region beta is zero", level=AssertionLevel.warning);
    beta := Modelica.Constants.small;
  else
    assert(false, "isobaricExpansionCoefficient error: phase is neither 1 nor 2");
  end if;
end isobaricExpansionCoefficient_cache;
