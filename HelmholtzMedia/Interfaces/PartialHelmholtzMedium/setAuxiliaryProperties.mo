within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setAuxiliaryProperties
  "intermediate results shared by the property functions of one state: Helmholtz derivatives (single phase) or saturation properties (two phase)"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 2 for two-phase, otherwise single phase";
  output AuxiliaryProperties aux;

algorithm
  if phase == 2 then
    aux.sat := setSat_T(T=T);
    aux.fl := EoS.setHelmholtzDerivs(d=aux.sat.liq.d, T=T, order=2);
    aux.fv := EoS.setHelmholtzDerivs(d=aux.sat.vap.d, T=T, order=2);
  else
    aux.f := EoS.setHelmholtzDerivs(d=d, T=T, order=2);
  end if;

  annotation (Documentation(info="<html>
<p>
Computes the intermediate results that the property functions of a state need: for a single-phase state the
derivatives of the Helmholtz energy at the state, for a two-phase state the saturation properties at the state's
temperature and the derivatives at the saturated liquid and vapour densities. Only the part for the phase of the
state is computed; the other fields are not set.
</p>
<p>
The property functions (specificHeatCapacityCp, velocityOfSound, density_derp_h, ...) are inlined one-line wrappers
that call their worker function (e.g. specificHeatCapacityCp_aux) with the result of this function. All of them call
setAuxiliaryProperties with the same arguments, so Dymola evaluates it only once per state and model evaluation.
The inputs are scalars on purpose: identical calls with the state record as argument are not merged by Dymola
when they come from different inlined functions.
</p>
</html>"));
end setAuxiliaryProperties;
