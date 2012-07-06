within HelmholtzMedia.Examples.Tests.MediaTestModels;
model ButaneTestModel "Test HelmholtzMedia.HelmholtzFluids.Butane"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=2.01, NumberOfIntervals=1000),
    __Dymola_experimentSetupOutput);
end ButaneTestModel;
