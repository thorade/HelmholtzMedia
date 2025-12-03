within HelmholtzMedia.Examples.MediaTestModels;
model IsopentaneTestModel "Test HelmholtzMedia.HelmholtzFluids.Isopentane"
  extends Modelica.Icons.Example;
  extends HelmholtzMedia.Examples.MediaTestModels.PartialTestModel(redeclare package Medium =
        HelmholtzMedia.HelmholtzFluids.Isopentane);

  annotation (experiment(StopTime=11));

end IsopentaneTestModel;
