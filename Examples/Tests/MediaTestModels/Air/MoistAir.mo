within HelmholtzMedia.Examples.Tests.MediaTestModels.Air;
model MoistAir "Test Modelica.Media.Air.MoistAir"
    extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = Modelica.Media.Air.MoistAir);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end MoistAir;
