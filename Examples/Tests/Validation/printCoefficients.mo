within HelmholtzMedia.Examples.Tests.Validation;
model printCoefficients "pretty printing of EoS coefficients"
  package medium = HelmholtzFluids.Ethanol;

algorithm
    // print to printlog.txt or to Simulation Log -> Simulation
    // if printing fails with error "Room to allocate string"
    // go to $Dymola$/source/matrixop.h and increase the size of simplestring by a factor of 10

    Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|", "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealLog), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealPower), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealEinstein), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualPoly), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualBwr), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualGauss), "printlog.txt");

    Modelica.Utilities.Streams.print("===============================================================================", "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.pressureSaturation), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.densityLiquid), "printlog.txt");
    Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.densityVapor), "printlog.txt");

  annotation (experiment(NumberOfIntervals=1), __Dymola_experimentSetupOutput);
end printCoefficients;
