within HelmholtzMedia.Examples.MediaTestModels;
model EthanolTestModel "Test HelmholtzMedia.HelmholtzFluids.Ethanol"
  extends Modelica.Icons.Example;
  extends HelmholtzMedia.Examples.MediaTestModels.PartialTestModel(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ethanol,
    volume(use_p_start=true, use_T_start=false),
    fixedMassFlowRate(use_T_ambient=false),
    ambient(use_T_ambient=false));

  annotation (experiment(StopTime=11));

end EthanolTestModel;
