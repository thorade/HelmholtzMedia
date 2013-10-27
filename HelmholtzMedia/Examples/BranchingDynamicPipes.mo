within HelmholtzMedia.Examples;
model BranchingDynamicPipes
  extends Modelica.Fluid.Examples.BranchingDynamicPipes(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane);
end BranchingDynamicPipes;
