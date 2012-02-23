within HelmholtzFluids.Examples.Tests.MediaTestModels.IdealGases;
model Air "Test single gas Modelica.Media.IdealGases.SingleGases.Air"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.Air.DryAirNasa);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end Air;
