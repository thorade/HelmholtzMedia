within HelmholtzMedia.Examples.Tests.MediaTestModels.Water;
model WaterIF97OnePhase_ph "Test Modelica.Media.Water.WaterIF97OnePhase_ph"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.Water.WaterIF97OnePhase_ph,
    fixedMassFlowRate(use_T_ambient=false, h_ambient=363755),
    ambient(use_T_ambient=false, h_ambient=112570));
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end WaterIF97OnePhase_ph;
