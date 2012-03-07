within HelmholtzMedia.Examples.Tests.MediaTestModels.HelmholtzFluids;
model R134a "Test Modelica.Media.Water.IdealSteam"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = Modelica.Media.Water.IdealSteam);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end R134a;
