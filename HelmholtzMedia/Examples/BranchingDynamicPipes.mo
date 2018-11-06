within HelmholtzMedia.Examples;
model BranchingDynamicPipes
  extends Modelica.Fluid.Examples.BranchingDynamicPipes(
    redeclare package Medium =
        HelmholtzMedia.HelmholtzFluids.Carbondioxide,
    system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial));
  annotation (experiment(
      StopTime=2,
      Tolerance=1e-04));
end BranchingDynamicPipes;
