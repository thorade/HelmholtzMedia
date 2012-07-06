within HelmholtzMedia.Examples.Tests.MediaTestModels;
model IsobutaneTestModel "Test HelmholtzMedia.HelmholtzFluids.Isobutane"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Isobutane);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end IsobutaneTestModel;
