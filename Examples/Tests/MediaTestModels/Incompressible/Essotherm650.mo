within HelmholtzFluids.Examples.Tests.MediaTestModels.Incompressible;
model Essotherm650 "Test Modelica.Media.Incompressible.Examples.Essotherm65"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.Incompressible.Examples.Essotherm650);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end Essotherm650;
