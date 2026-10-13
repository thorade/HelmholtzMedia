within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
package Auxiliary "Internal functions: shared intermediate results, worker functions of the property functions, and state/derivative functions"
  extends Modelica.Icons.InternalPackage;

  annotation (Documentation(info="<html>
<p>
Internal functions of PartialHelmholtzMedium; they are not meant to be called directly in models.
</p>
<ul>
<li><b>AuxiliaryProperties</b> and <b>setAuxiliaryProperties</b>: intermediate results shared by the property functions of one state
(Helmholtz derivatives in single phase, saturation properties in two phase).</li>
<li><b>*_aux</b>: worker functions of the property functions of a state (e.g. specificHeatCapacityCp_aux).
The public property function is an inlined one-line wrapper that calls the worker with setAuxiliaryProperties(state.d, state.T, state.phase).
All wrappers call setAuxiliaryProperties with the same scalar arguments, so the tool evaluates it only once per state.</li>
<li><b>*_state</b>: property functions of two independent variables (e.g. density_ph_state) that receive the thermodynamic state as an
additional argument and carry the derivative annotation; the public function (e.g. density_ph) keeps the inverse annotation.</li>
<li><b>*_der</b>: time derivatives used by the derivative annotations of the *_state functions.</li>
</ul>
</html>"));
end Auxiliary;
