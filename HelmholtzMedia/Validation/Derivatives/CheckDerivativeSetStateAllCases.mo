within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeSetStateAllCases "CheckDerivativeSetState for all end states"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
  CheckDerivativeSetState singlePhase(redeclare package Medium = Medium, endStates=EndStates.SinglePhase);
  CheckDerivativeSetState twoPhase(redeclare package Medium = Medium, endStates=EndStates.TwoPhase);
  CheckDerivativeSetState phaseChange(redeclare package Medium = Medium, endStates=EndStates.PhaseChange);
  CheckDerivativeSetState liquidNearSaturation(redeclare package Medium = Medium, endStates=EndStates.LiquidNearSaturation);
  CheckDerivativeSetState vapourNearSaturation(redeclare package Medium = Medium, endStates=EndStates.VapourNearSaturation);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeSetStateAllCases;
