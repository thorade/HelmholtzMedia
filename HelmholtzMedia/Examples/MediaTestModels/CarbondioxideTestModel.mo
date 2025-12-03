within HelmholtzMedia.Examples.MediaTestModels;
model CarbondioxideTestModel "Test HelmholtzMedia.HelmholtzFluids.Carbondioxide"
  extends Modelica.Icons.Example;
  extends HelmholtzMedia.Examples.MediaTestModels.PartialTestModel(redeclare package Medium =
        HelmholtzMedia.HelmholtzFluids.Carbondioxide);

  annotation (experiment(StopTime=11));

end CarbondioxideTestModel;
