within HelmholtzMedia.Examples.Tests.MediaTestModels.Water;
model ConstantPropertyLiquidWater
  "Test Modelica.Media.Water.ConstantPropertyLiquidWater"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.Water.ConstantPropertyLiquidWater);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end ConstantPropertyLiquidWater;
