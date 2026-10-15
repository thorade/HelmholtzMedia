within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsPentane "CheckDerivativeAnnotationsAllCases for Pentane"
  extends CheckDerivativeAnnotationsAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Pentane);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsPentane;
