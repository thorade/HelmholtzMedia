within HelmholtzMedia.Examples.MediaTestModels;
model ButaneTestModel_dT_component_ph "Test HelmholtzMedia.HelmholtzFluids.Butane"
  extends Modelica.Icons.Example;
  extends HelmholtzMedia.Examples.MediaTestModels.PartialTestModel(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane (inputChoice=Medium.InputChoice.dT),
    ambient(use_p_ambient=true, use_T_ambient=false),
    volume(use_T_start=false, medium(preferredMediumStates=true, componentInputChoice=Medium.InputChoice.ph)));

  annotation (experiment(StopTime=1.01), Diagram(graphics));

end ButaneTestModel_dT_component_ph;
