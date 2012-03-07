within HelmholtzMedia.Examples.Tests.MediaTestModels.Water;
model WaterIF97_pT "Test Modelica.Media.Water.WaterIF97_pT"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = Modelica.Media.Water.WaterIF97_pT);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end WaterIF97_pT;
