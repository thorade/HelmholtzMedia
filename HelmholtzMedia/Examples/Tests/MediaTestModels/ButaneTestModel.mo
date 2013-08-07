within HelmholtzMedia.Examples.Tests.MediaTestModels;
model ButaneTestModel "Test HelmholtzMedia.HelmholtzFluids.Butane"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane);

  annotation (experiment(StopTime=1.01));

end ButaneTestModel;
