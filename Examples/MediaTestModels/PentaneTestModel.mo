within HelmholtzMedia.Examples.MediaTestModels;
model PentaneTestModel "Test HelmholtzMedia.HelmholtzFluids.Pentane"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Tests.Components.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Pentane,
    fixedMassFlowRate(use_T_ambient=false),
    volume(use_T_start=false),
    ambient(use_T_ambient=false),
    h_start=200000,
    p_start=1000000);

  annotation (experiment(StopTime=1000));

end PentaneTestModel;
