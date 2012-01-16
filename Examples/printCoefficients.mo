within HelmholtzFluids.Examples;
model printCoefficients
  package medium = HelmholtzFluids.Butane;

algorithm
    // print Coefficients to Simulation Log -> Simulation
    Modelica.Utilities.Streams.print("=================================================");
    Modelica.Utilities.Streams.print(Modelica.Math.Vectors.toString(medium.helmholtzCoefficients.n_ideal));
    Modelica.Utilities.Streams.print(Modelica.Math.Vectors.toString(medium.helmholtzCoefficients.n_residual));

  annotation (experiment(NumberOfIntervals=1), __Dymola_experimentSetupOutput);
end printCoefficients;
