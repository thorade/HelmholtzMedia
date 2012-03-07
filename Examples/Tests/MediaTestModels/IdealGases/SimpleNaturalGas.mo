within HelmholtzMedia.Examples.Tests.MediaTestModels.IdealGases;
model SimpleNaturalGas
  "Test mixture gas Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGas"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGas);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end SimpleNaturalGas;
