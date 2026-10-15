within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
package Internal "Internal implementation: equation of state, ancillary equations, transport models, shared cache and worker functions"
  extends Modelica.Icons.InternalPackage;

  annotation (Documentation(info="<html>
<p>
Internal implementation of PartialHelmholtzMedium. The functions promised by the medium interface
(Modelica.Media.Interfaces.PartialTwoPhaseMedium) and the additional property functions of HelmholtzMedia
are located directly in the medium package; the functions here are not meant to be called in models.
</p>
<ul>
<li><b>EoS</b>: Helmholtz energy equation of state (derivatives of the Helmholtz energy and the thermodynamic properties computed from them).</li>
<li><b>Ancillary</b>: ancillary equations (saturation pressure and densities, melting curve, start values).</li>
<li><b>Transport</b>: transport property models (viscosity, thermal conductivity, surface tension).</li>
<li><b>HelmholtzCache</b> and <b>setHelmholtzCache</b>: intermediate results shared by the property functions of one state
(Helmholtz derivatives in single phase, saturation properties in two phase).</li>
<li><b>*_cache</b>: worker functions of the property functions (e.g. specificHeatCapacityCp_cache).
The public property function is an inlined one-line wrapper that calls the worker with setHelmholtzCache(state.d, state.T, state.phase).
All wrappers call setHelmholtzCache with the same scalar arguments, so the tool evaluates it only once per state.</li>
<li><b>*_state</b>: property functions of two independent variables (e.g. density_ph_state) that receive the thermodynamic state as an
additional argument and carry the derivative annotation; the public function (e.g. density_ph) keeps the inverse annotation.</li>
<li><b>*_state_der</b>: time derivatives of the *_state functions (their derivative annotation).</li>
<li><b>*_der</b>: time derivatives of the public functions of two variables (their derivative annotation, used e.g. for
analytic Jacobians); inlined, they call the *_state_der function with the same setState call as the public function.</li>
<li><b>*_cache_der</b>: time derivatives of the worker functions (derivative annotation with noDerivative=cache). The workers
take the scalars d, T and phase instead of the state record, so that the tool can differentiate them; the derivatives need
the Helmholtz derivatives up to third order, which are computed in the *_cache_der function.
Available for all workers; the saturation workers (dBubbleDensity_dPressure_cache etc.) take sat.Tsat as additional scalar
input, since along the saturation line they depend on T only.</li>
<li><b>setState_*X_der</b>, <b>setSat_p_der</b>, <b>setSat_T_der</b>: derivative annotations of setState_phX, setState_pTX,
setState_dTX, setState_psX, setSat_p and setSat_T, so that functions of a state or saturation state computed by them
(as in Modelica.Fluid components) can be differentiated, e.g. for analytic Jacobians. The derivative has the same record
type as the original output (phase=0). (The Modelica Language Specification 3.8-dev, section 12.7.1, describes a derivative
record with the Real components only; Dymola 2027x expects the same record type.)
Derivatives along the saturation line use Clausius-Clapeyron and its derivative, see Thorade and Saadat (2013).</li>
<li><b>setSat_T_Newton</b>: Newton iteration of setSat_T (start values from Ancillary.setSat_T_start).</li>
</ul>
</html>"));
end Internal;
