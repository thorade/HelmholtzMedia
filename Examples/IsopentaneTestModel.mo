within HelmholtzMedia.Examples;
model IsopentaneTestModel "Test HelmholtzMedia.HelmholtzFluids.Isopentane"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Isopentane);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end IsopentaneTestModel;
