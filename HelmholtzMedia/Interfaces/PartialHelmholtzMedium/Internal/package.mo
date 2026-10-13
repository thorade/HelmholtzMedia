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
<li><b>*_der</b>: time derivatives used by the derivative annotations of the *_state functions.</li>
<li><b>setSat_T_Newton</b>: Newton iteration of setSat_T (start values from Ancillary.setSat_T_start).</li>
</ul>
</html>"));
end Internal;
