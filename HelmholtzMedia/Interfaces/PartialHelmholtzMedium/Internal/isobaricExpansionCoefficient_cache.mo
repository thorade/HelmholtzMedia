within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function isobaricExpansionCoefficient_cache "returns 1/v*(dv/dT)@p=const (computed from the shared auxiliary properties, see isobaricExpansionCoefficient)"
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  output IsobaricExpansionCoefficient beta;
algorithm
  if phase==1 then
    beta := 1.0/d*EoS.dpTd(cache.f)/EoS.dpdT(cache.f);
  elseif phase==2 then
    assert(false, "isobaricExpansionCoefficient warning: in the two-phase region beta is zero", level=AssertionLevel.warning);
    beta := Modelica.Constants.small;
  else
    assert(false, "isobaricExpansionCoefficient error: phase is neither 1 nor 2");
  end if;
  annotation (derivative(noDerivative=cache)=isobaricExpansionCoefficient_cache_der);
end isobaricExpansionCoefficient_cache;
