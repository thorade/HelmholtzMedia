within HelmholtzMedia.Examples.Tests.MediaTestModels;
model R134aTestModel "Test HelmholtzMedia.HelmholtzFluids.R134a"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a);

  annotation (experiment(StopTime=11));

end R134aTestModel;
