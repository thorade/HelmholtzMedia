within HelmholtzFluids.Examples.Tests.MediaTestModels.IdealGases;
model Nitrogen "Test single gas Modelica.Media.IdealGases.SingleGases.N2"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.IdealGases.SingleGases.N2);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end Nitrogen;
