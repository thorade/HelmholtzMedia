within HelmholtzMedia.HelmholtzFluids;
package Isopentane "Isopentane"
  extends Interfaces.PartialHelmholtzMedium(
    fluidConstants={fluidConstantsIsopentane},
    helmholtzCoefficients=helmholtzCoefficientsIsopentane,
    thermalConductivityCoefficients=thermalConductivityCoefficientsIsopentane,
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsIsopentane,
    surfaceTensionCoefficients=surfaceTensionCoefficientsIsopentane,
    ancillaryCoefficients=ancillaryCoefficientsIsopentane,
    fluidLimits=fluidLimitsIsopentane,
    Density(min=fluidLimitsIsopentane.DMIN, max=fluidLimitsIsopentane.DMAX, start=fluidConstantsIsopentane.molarMass/fluidConstantsIsopentane.criticalMolarVolume),
    Temperature(min=fluidLimitsIsopentane.TMIN, max=fluidLimitsIsopentane.TMAX, start=298.15),
    AbsolutePressure(min=0, max=35e6, start=101325),
    SpecificEnthalpy(min=fluidLimitsIsopentane.HMIN, max=fluidLimitsIsopentane.HMAX, start=(fluidLimitsIsopentane.HMIN+fluidLimitsIsopentane.HMAX)/2),
    SpecificEntropy(min=fluidLimitsIsopentane.SMIN, max=fluidLimitsIsopentane.SMAX, start=(fluidLimitsIsopentane.SMIN+fluidLimitsIsopentane.SMAX)/2));

  final constant FluidConstants
  fluidConstantsIsopentane(
    chemicalFormula="C5H12",
    structureFormula="",
    casRegistryNumber="78-78-4",
    iupacName="",
    molarMass=0.07214878,
    hasCriticalData=true,
       criticalTemperature=460.35,
       criticalPressure=3378000,
       criticalMolarVolume=1/3271,
       HCRIT0=503558.750553544,
       SCRIT0=1282.74524346652,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.11,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.2274,
    triplePointTemperature=112.65,
    triplePointPressure=0.000089527,
    normalBoilingPoint=300.98,
    meltingPoint=112.662) "Fluid Constants";

  final constant FluidLimits
  fluidLimitsIsopentane(
    TMIN=fluidConstantsIsopentane.triplePointTemperature,
    TMAX=500,
    DMIN=Modelica.Constants.small,
    DMAX=960,
    PMIN=Modelica.Constants.small,
    PMAX=1000e6,
    HMIN=-725e3,
    HMAX=+700e3,
    SMIN=-3036,
    SMAX=9283) "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsIsopentane(
    idealLog=[
      +3.0000000000,    1],
    idealPower=[
      +2.5822330405,    0;
       1.1609103419,    1],
    idealEinstein=[
      +7.4056000000,   -0.9601390247;
       9.5772000000,   -2.4090366026;
      15.7650000000,   -4.4944064299;
      12.1190000000,   -9.1082871728],
    residualPoly=[
      +1.0963,          0.25,    1.0;
      -3.0402,          1.125,   1.0;
       1.0317,          1.5,     1.0;
      -0.15410,         1.375,   2.0;
       0.11535,         0.25,    3.0;
       0.00029809,      0.875,   7.0],
    residualBwr=[
      +0.39571,         0.625,   2.0,   1;
      -0.045881,        1.75,    5.0,   1;
      -0.35804,         3.625,   1.0,   2;
      -0.10107,         3.625,   4.0,   2;
      -0.035484,       14.5,     3.0,   3;
       0.018156,       12.0,     4.0,   3],
   residualGauss=fill(0.0, 0, 9)) "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsIsopentane(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0 = 341.06,
    reducingThermalConductivity_0 = 1e-3,
    lambda_0_num_coeffs=[
       1.35558587,              0.0;
      -0.152666315743857,      -1.0;
       1.,                    -96.0],
    reducingTemperature_background=460.51,
    reducingMolarVolume_background=1/3240,
    reducingThermalConductivity_background=1e-3,
    lambda_b_coeffs=[
      18.608933103800,  0.0,    1.0,   0.0;
      -5.836570612990,  0.0,    3.0,   0.0;
       3.489871005290,  0.0,    4.0,   0.0;
       0.704467355508, -1.0,    4.0,   0.0;
      -0.206501417728,  0.0,    5.0,   0.0;
      -0.223070394020, -1.0,    5.0,   0.0],
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.9316E-9,
    T_ref=690.525) "Coefficients for the thermal conductivity";

  final constant Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsIsopentane(
    dynamicViscosityModel=DynamicViscosityModel.VS2,
    collisionIntegralModel=CollisionIntegralModel.CI0,
    sigma=0.56232,
    epsilon_kappa=341.06,
    CET=[
      0.2267237, 0.5],
    b=[
      0.0,    0.0;
      0.0,    0.0;
      0.0,    0.0;
      100.0,  0.0],
    c=[
      -4.57981980159405;
      -3393.52438560000;
       9.38066543240000;
       33641.3512000000;
       0.15624235969000;
       122.900175430000;
      -20914.7951660000;
       3.24]) "Coefficients for the dynamic viscosity";

  final constant Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsIsopentane(
    coeffs=[
       0.05106,     1.21]) "Coefficients for the surface tension";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsIsopentane(
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.72392E+01,       1.0;
       0.22635E+01,       1.5;
      -0.18237E+01,       2.02;
      -0.29997E+01,       4.24;
      -0.27752E+01,      16.1],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       0.18367E+02,       1.21;
      -0.30283E+02,       1.41;
       0.13557E+02,       1.65;
      -0.90533E+00,       0.09;
       0.20927E+01,       0.164],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -0.38825E+02,       0.565;
       0.79040E+02,       0.66;
      -0.48791E+02,       0.77;
      -0.21603E+02,       3.25;
      -0.57218E+02,       7.3;
      -0.15164E+03,      16.6],
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=112.65,
    p_reducing=0.84e-4,
    pressureMelting1=[
      -7127700000000.0,       0;
       7127700000001.0,       1.563],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2))
  "Coefficients for the ancillary equations (PS5, DL2, DV6, ML1)";


  annotation (Documentation(info="<html>
These are the coefficients for Isopentane.

<dl>
<dt> Lemmon, E. W. and Span, R</dt>
<dd> <b>Short Fundamental Equations of State for 20 Industrial Fluids</b><br />
     Journal of Chemical &amp; Engineering Data 51 (3) , 785-850 (2006)<br />
     DOI: <a href=\"http://dx.doi.org/10.1021/je050186n\">10.1021/je050186n</a>
</dd>
<dt> Reeves, L.E. and Scott, G.J. and Babb, S.E.Jr.</dt>
<dd> <b>Melting curves of pressure-transmitting fluids</b><br />
     Journal of Chemical Physics 40 (12) , 3662-3666 (1964)<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.1725068\">10.1063/1.1725068</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Isopentane;
