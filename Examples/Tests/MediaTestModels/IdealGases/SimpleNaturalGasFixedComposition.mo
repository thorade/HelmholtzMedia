within HelmholtzFluids.Examples.Tests.MediaTestModels.IdealGases;
model SimpleNaturalGasFixedComposition
  "Test mixture gas Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGas"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium =
        Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGasFixedComposition);
  annotation (experiment(StopTime=1.01));
end SimpleNaturalGasFixedComposition;
