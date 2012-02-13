within HelmholtzFluids.Examples;
model printCoefficients "pretty printing of EoS coefficients"
  package medium = HelmholtzFluids.Butane;

algorithm
    // print to Simulation Log -> Simulation or print to printlog.txt
    // if printing fails with error "Room to allocate string"
    // go to $Dymola$/source/matrixop.h and increase the size of simplestring by a factor of 10

    Modelica.Utilities.Streams.print("===============================================================================", "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.ideal), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualPoly), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualBwr), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualGauss), "printlog.txt");

    Modelica.Utilities.Streams.print("===============================================================================", "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.pressureSaturation), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.densityLiquid), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.densityVapor), "printlog.txt");

  annotation (experiment(NumberOfIntervals=1), __Dymola_experimentSetupOutput);
end printCoefficients;
