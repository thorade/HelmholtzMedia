within HelmholtzMedia.Examples;
model BranchingDynamicPipes
  extends Modelica.Fluid.Examples.BranchingDynamicPipes(
    redeclare package Medium =
        HelmholtzMedia.HelmholtzFluids.Propane,
    system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial));

  annotation (experiment(
      StopTime=1,
      Tolerance=1e-05));
end BranchingDynamicPipes;
