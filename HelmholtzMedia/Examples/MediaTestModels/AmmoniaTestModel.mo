within HelmholtzMedia.Examples.MediaTestModels;
model AmmoniaTestModel "Test HelmholtzMedia.HelmholtzFluids.Ammonia"
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Utilities.PartialTestModel(
     redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia,
      fixedMassFlowRate(use_T_ambient=false));

  annotation (experiment(StopTime=11));

end AmmoniaTestModel;
