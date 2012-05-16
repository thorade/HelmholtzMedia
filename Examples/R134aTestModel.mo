within HelmholtzMedia.Examples;
model R134aTestModel "Test HelmholtzMedia.HelmholtzFluids.R134a"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a(independentVariables=IndependentVariables.ph));
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=1.01));
end R134aTestModel;
