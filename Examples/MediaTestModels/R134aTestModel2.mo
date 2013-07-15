within HelmholtzMedia.Examples.MediaTestModels;
model R134aTestModel2 "Test HelmholtzMedia.HelmholtzFluids.R134a"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel2(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a);

  annotation (experiment(StopTime=11));

end R134aTestModel2;
