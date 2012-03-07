within HelmholtzMedia.Examples.Tests.MediaTestModels.Air;
model DryAirNasa "Test Modelica.Media.Air.DryAirNasa"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = Modelica.Media.Air.DryAirNasa);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end DryAirNasa;
