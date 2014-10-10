within HelmholtzMedia.Examples.Validation;
model printCoefficients "pretty printing of EoS coefficients"
  package Medium = HelmholtzFluids.Helium;

protected
  String fileName = "printCoefficients.txt";
  final constant Real kilo = 1e3;

algorithm
  when terminal() then
  // if printing fails with error "Room to allocate string"
  // go to $Dymola$/source/matrixop.h and increase the size of simplestring, e.g. by a factor of 10

  // remove old file
  Modelica.Utilities.Files.remove(fileName);

  Modelica.Utilities.Streams.print(Medium.mediumName, fileName);
  Modelica.Utilities.Streams.print(Medium.fluidConstants[1].casRegistryNumber, fileName);
  Modelica.Utilities.Streams.print(Medium.fluidConstants[1].iupacName, fileName);
  Modelica.Utilities.Streams.print(Medium.fluidConstants[1].structureFormula + " ----> " + Medium.fluidConstants[1].chemicalFormula, fileName);
  Modelica.Utilities.Streams.print("", fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].molarMass*kilo, significantDigits=9), fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].triplePointTemperature), fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].normalBoilingPoint), fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].criticalTemperature), fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].criticalPressure/kilo), fileName);
  Modelica.Utilities.Streams.print(String(1.0/Medium.fluidConstants[1].criticalMolarVolume/kilo, significantDigits=10), fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].acentricFactor), fileName);
  Modelica.Utilities.Streams.print(String(Medium.fluidConstants[1].dipoleMoment), fileName);

  Modelica.Utilities.Streams.print("\n =============================================================================== \n #FEQ", fileName); // 80 characters
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.residualPoly), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.residualBwr), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.residualGauss), fileName);

  Modelica.Utilities.Streams.print("\n =============================================================================== \n #AUX: CPP or PH0", fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.idealLog), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.idealPower), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.idealEinstein), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.idealCosh), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.helmholtzCoefficients.idealSinh), fileName);

  Modelica.Utilities.Streams.print("\n =============================================================================== \n #STN", fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.surfaceTensionCoefficients.coeffs), fileName);
  Modelica.Utilities.Streams.print("#MLT", fileName);
  Modelica.Utilities.Streams.print(String(Medium.ancillaryCoefficients.pressureMeltingModel), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.ancillaryCoefficients.pressureMelting1), fileName);
//Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.ancillaryCoefficients.pressureMelting2), fileName);
//Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.ancillaryCoefficients.pressureMelting3), fileName);
  Modelica.Utilities.Streams.print("#PS", fileName);
  Modelica.Utilities.Streams.print(String(Medium.ancillaryCoefficients.pressureSaturationModel), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.ancillaryCoefficients.pressureSaturation), fileName);
  Modelica.Utilities.Streams.print("#DL", fileName);
  Modelica.Utilities.Streams.print(String(Medium.ancillaryCoefficients.densityLiquidModel), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.ancillaryCoefficients.densityLiquid), fileName);
  Modelica.Utilities.Streams.print("#DV", fileName);
  Modelica.Utilities.Streams.print(String(Medium.ancillaryCoefficients.densityVaporModel), fileName);
  Modelica.Utilities.Streams.print(Modelica.Math.Matrices.toString(Medium.ancillaryCoefficients.densityVapor), fileName);
  Modelica.Utilities.Streams.print("===============================================================================", fileName);
  end when;
end printCoefficients;
