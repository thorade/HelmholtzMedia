within HelmholtzMedia.HelmholtzFluids;
package Pentane "Pentane"
  extends Interfaces.PartialHelmholtzMedium(
    fluidConstants={fluidConstantsPentane},
    helmholtzCoefficients=helmholtzCoefficientsPentane,
    thermalConductivityCoefficients=thermalConductivityCoefficientsPentane,
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsPentane,
    surfaceTensionCoefficients=surfaceTensionCoefficientsPentane,
    ancillaryCoefficients=ancillaryCoefficientsPentane,
    fluidLimits=fluidLimitsPentane,
    Density(min=fluidLimitsPentane.DMIN, max=fluidLimitsPentane.DMAX, start=fluidConstantsPentane.molarMass/fluidConstantsPentane.criticalMolarVolume),
    Temperature(min=fluidLimitsPentane.TMIN, max=fluidLimitsPentane.TMAX, start=298.15),
    AbsolutePressure(min=Modelica.Constants.small, max=100e6, start=101325),
    SpecificEnthalpy(min=fluidLimitsPentane.HMIN, max=fluidLimitsPentane.HMAX, start=(fluidLimitsPentane.HMIN+fluidLimitsPentane.HMAX)/2),
    SpecificEntropy(min=fluidLimitsPentane.SMIN, max=fluidLimitsPentane.SMAX, start=(fluidLimitsPentane.SMIN+fluidLimitsPentane.SMAX)/2));

  final constant Interfaces.PartialHelmholtzMedium.FluidConstants
  fluidConstantsPentane(
    chemicalFormula="C5H12",
    structureFormula="",
    casRegistryNumber="109-66-0",
    iupacName="",
    molarMass=0.07214878,
    hasCriticalData=true,
       criticalTemperature=469.7,
       criticalPressure=3370000,
       criticalMolarVolume=1/3215.5776,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.07,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.251,
    triplePointTemperature=143.47,
    triplePointPressure=0.000076322,
    normalBoilingPoint=309.21,
    meltingPoint=143.47) "Fluid Constants";

  final constant Interfaces.PartialHelmholtzMedium.FluidLimits
  fluidLimitsPentane(
    TMIN=fluidConstantsPentane.triplePointTemperature,
    TMAX=600,
    DMIN=Modelica.Constants.small,
    DMAX=808.066,
    PMIN=Modelica.Constants.small,
    PMAX=100e6,
    HMIN=-200e3,
    HMAX=+1500e3,
    SMIN=-1000,
    SMAX=5000) "Helmholtz EoS Limits";

  final constant Interfaces.PartialHelmholtzMedium.EoS.HelmholtzCoefficients
  helmholtzCoefficientsPentane(
    idealLog=[
      3,    1],
    idealPower=[
      7.840639316,     0;
      -84.68310031,    1],
    idealEinstein=fill(0.0, 0, 2),
    idealCosh=[
      21.836,   1.789520971],
    idealSinh=[
      8.95043,    0.380391739;
      33.4032,    3.777411113],
    residualPoly=[
      1.0968643,        0.25,    1.0,   0;
     -2.9988888,        1.125,   1.0,   0;
      0.99516887,       1.5,     1.0,   0;
     -0.16170709,       1.375,   2.0,   0;
      0.1133446,        0.25,    3.0,   0;
      0.000267606,      0.875,   7.0,   0],
    residualBwr=[
      0.40979882,         0.625,   2.0,   1;
     -0.040876423,        1.75,    5.0,   1;
     -0.38169482,         3.625,   1.0,   2;
     -0.10931957,         3.625,   4.0,   2;
     -0.032073223,       14.5,     3.0,   3;
      0.016877016,       12.0,     4.0,   3],
   residualGauss=fill(0.0, 0, 9)) "Coefficients of the Helmholtz EoS";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsPentane(
    reducingTemperature_0 = 341.1,
    reducingThermalConductivity_0 = 1e-3,
    lambda_0_num_coeffs=[
       1.35558587,             0.0;
      -0.15569137,            -1.0;
       1.,                    -96.0],
    reducingTemperature_residual=469.69,
    reducingMolarVolume_residual=1/3215,
    reducingThermalConductivity_residual=1e-3,
    lambda_r_coeffs=[
      18.608933103800,  0.0,    1.0,   0.0;
      -5.836570612990,  0.0,    3.0,   0.0;
       3.489871005290,  0.0,    4.0,   0.0;
       0.704467355508, -1.0,    4.0,   0.0;
      -0.206501417728,  0.0,    5.0,   0.0;
      -0.223070394020, -1.0,    5.0,   0.0],
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.9345E-9,
    T_ref=704.55) "Coefficients for the thermal conductivity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsPentane(
    dynamicViscosityModel=Interfaces.PartialHelmholtzMedium.Types.DynamicViscosityModel.VS2,
    collisionIntegralModel=Interfaces.PartialHelmholtzMedium.Types.CollisionIntegralModel.CI0,
    sigma=0.5784,
    epsilon_kappa=341.10,
    CET=[
      0.226720214, 0.5],
    b=[
      0.0,    0.0;
      0.0,    0.0;
      0.0,    0.0;
      100.0,  0.0],
    c=[
      -13.47938293;
       1176.62751650;
       14.2278439927;
      -21951.0293411;
       0.03766867689;
       70.1529173825;
       21435.7720323;
       3.215]) "Coefficients for the dynamic viscosity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsPentane(
    coeffs=[
     0.0562267,   1.25;
    -0.0037496,   2.25;
    -0.0029861,   3.25]) "Coefficients for the surface tension";

  final constant
  Interfaces.PartialHelmholtzMedium.Ancillary.AncillaryCoefficients
  ancillaryCoefficientsPentane(
    pressureSaturationModel=Interfaces.PartialHelmholtzMedium.Types.PressureSaturationModel.PS5,
    pressureSaturation=[
     -0.73918E+01,   1.0;
      0.31102E+01,   1.5;
     -0.22415E+01,   1.74;
     -0.31585E+01,   3.75;
     -0.90451E+00,   8.0],
    densityLiquidModel=Interfaces.PartialHelmholtzMedium.Types.DensityLiquidModel.DL1,
    densityLiquid=[
      0.10178E+01,   0.27;
      0.42703E+00,   0.44;
      0.11334E+01,   0.6;
      0.41518E+00,   4.0;
     -0.47950E-01,   5.0],
    densityVaporModel=Interfaces.PartialHelmholtzMedium.Types.DensityVaporModel.DV3,
    densityVapor=[
     -0.29389E+01,   0.4;
     -0.62784E+01,   1.18;
     -0.19941E+02,   3.2;
     -0.16709E+02,   6.6;
     -0.36543E+02,   7.0;
     -0.12799E+03,  15.0])
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for Pentane.

<dl>
<dt> Span, R. & Wagner, W. </dt>
<dd> <b>Equations of State for Technical Applications. II. Results for Nonpolar Fluids</b><br>
     International Journal of Thermophysics, Vol. 24, No. 1, (2003), pp. 41-109 <br>
     DOI: <a href=\"http://dx.doi.org/10.1023/A:1022310214958\">10.1023/A:1022310214958</a>
</dd>
<dt> Jaeschke, M. & Schley, P. </dt>
<dd> <b>Ideal-gas thermodynamic properties for natural-gas applications</b><br>
     International Journal of Thermophysics, Vol. 16, No. 6 (1995), pp. 1381-1392<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF02083547\">10.1007/BF02083547</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Pentane;
