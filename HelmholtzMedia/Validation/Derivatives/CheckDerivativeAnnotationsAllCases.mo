within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsAllCases
  "CheckDerivativeAnnotations for all end states (single phase, two phase, phase change, close to saturation)"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
  CheckDerivativeAnnotations singlePhase(redeclare package Medium = Medium, endStates=EndStates.SinglePhase);
  CheckDerivativeAnnotations twoPhase(redeclare package Medium = Medium, endStates=EndStates.TwoPhase);
  CheckDerivativeAnnotations phaseChange(redeclare package Medium = Medium, endStates=EndStates.PhaseChange);
  CheckDerivativeAnnotations liquidNearSaturation(redeclare package Medium = Medium, endStates=EndStates.LiquidNearSaturation);
  CheckDerivativeAnnotations vapourNearSaturation(redeclare package Medium = Medium, endStates=EndStates.VapourNearSaturation);
  annotation (experiment(StopTime=1, Tolerance=1e-10), Documentation(info="<html>
<p>
Runs <a href=\"modelica://HelmholtzMedia.Validation.Derivatives.CheckDerivativeAnnotations\">CheckDerivativeAnnotations</a>
for all <a href=\"modelica://HelmholtzMedia.Validation.Derivatives.EndStates\">EndStates</a> of one fluid.
At the end of the simulation the relative deviation of each function is printed and an assert checks it against
the tolerance. Variants for other fluids redeclare <code>Medium</code>.
</p>
</html>"));
end CheckDerivativeAnnotationsAllCases;
