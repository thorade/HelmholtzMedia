within HelmholtzMedia.Examples.MediaTestModels;
model PropyleneTestModel "Test HelmholtzMedia.HelmholtzFluids.Propylene"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Utilities.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Propylene,
      fixedMassFlowRate(use_T_ambient=false));

  annotation (experiment(StopTime=11));

end PropyleneTestModel;
