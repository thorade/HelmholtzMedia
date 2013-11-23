within HelmholtzMedia.Examples;
model BranchingDynamicPipes
  extends Modelica.Fluid.Examples.BranchingDynamicPipes(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane, system(
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial));
end BranchingDynamicPipes;
