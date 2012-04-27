within HelmholtzMedia.Examples;
model R134aTestModel "Test Modelica.Media.Water.IdealSteam"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end R134aTestModel;
