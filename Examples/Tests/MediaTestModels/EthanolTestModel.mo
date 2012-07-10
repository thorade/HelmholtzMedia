within HelmholtzMedia.Examples.Tests.MediaTestModels;
model EthanolTestModel "Test HelmholtzMedia.HelmholtzFluids.Ethanol"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ethanol);

  annotation (experiment(StopTime=11));

end EthanolTestModel;
