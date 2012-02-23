within HelmholtzFluids.Examples.Tests.MediaTestModels.LinearFluid;
model LinearWater_pT "Test Modelica.Media.Incompressible.Examples.Essotherm65"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.CompressibleLiquids.LinearWater_pT_Ambient);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end LinearWater_pT;
