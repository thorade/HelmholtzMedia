within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsTwoPhase "as CheckDerivativeAnnotations, with two-phase end states"
  extends CheckDerivativeAnnotations(twoPhase=true);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsTwoPhase;
