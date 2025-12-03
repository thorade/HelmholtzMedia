within HelmholtzMedia.Examples.MediaTestModels;
model HeliumTestModel "Test HelmholtzMedia.HelmholtzFluids.Helium"
  extends Modelica.Icons.Example;
  extends HelmholtzMedia.Examples.MediaTestModels.PartialTestModel(redeclare package Medium =
        HelmholtzMedia.HelmholtzFluids.Helium, ambient(use_p_ambient=true, use_T_ambient=false));

  annotation (experiment(StopTime=1.01), Diagram(graphics));

end HeliumTestModel;
