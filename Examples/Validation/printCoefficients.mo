within HelmholtzMedia.Examples.Validation;
model printCoefficients "pretty printing of EoS coefficients"
  package medium = HelmholtzFluids.Helium;

protected
  String fileName = "printlog.txt";

algorithm
  // if printing fails with error "Room to allocate string"
  // go to $Dymola$/source/matrixop.h and increase the size of simplestring by a factor of 10

  // remove old file
  Modelica.Utilities.Files.remove(fileName);

  Modelica.Utilities.Streams.print("====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|====|", fileName); // 80 characters
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealLog), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealPower), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealEinstein), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealCosh), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.idealSinh), fileName);

  Modelica.Utilities.Streams.print("===============================================================================", fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualPoly), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualBwr), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.helmholtzCoefficients.residualGauss), fileName);

  Modelica.Utilities.Streams.print("===============================================================================", fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.pressureSaturation), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.densityLiquid), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(medium.ancillaryCoefficients.densityVapor), fileName);

annotation (experiment(NumberOfIntervals=1));
end printCoefficients;
