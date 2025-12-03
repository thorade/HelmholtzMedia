within HelmholtzMedia.Examples.MediaTestModels;
model R134aTestModel "Test HelmholtzMedia.HelmholtzFluids.R134a"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Utilities.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a,
     ambient(use_T_ambient=false),
     fixedMassFlowRate(use_T_ambient=false),
     volume(use_T_start=false));

  annotation (experiment(StopTime=11));

end R134aTestModel;
