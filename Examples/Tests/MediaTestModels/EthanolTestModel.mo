within HelmholtzMedia.Examples.Tests.MediaTestModels;
model EthanolTestModel "Test HelmholtzMedia.HelmholtzFluids.Ethanol"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ethanol);
  annotation (Documentation(info="<html>

</html>"),
   experiment(StopTime=2.01, NumberOfIntervals=1000),
    __Dymola_experimentSetupOutput);
end EthanolTestModel;
