within HelmholtzFluids.Examples;
model printCoefficients
  package medium = HelmholtzFluids.Butane;

algorithm
    // print Coefficients to Simulation Log -> Simulation
    Modelica.Utilities.Streams.print("=================================================");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.ideal));
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualPoly));
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualBwr));
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualGauss));

  annotation (experiment(NumberOfIntervals=1), __Dymola_experimentSetupOutput);
end printCoefficients;
