within HelmholtzFluids.Examples.Tests.MediaTestModels.Air;
model SimpleAir "Test Modelica.Media.Air.SimpleAir"
  extends Modelica.Icons.Example;
    extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = Modelica.Media.Air.SimpleAir);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end SimpleAir;
