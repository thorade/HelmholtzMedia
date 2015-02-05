within HelmholtzMedia.HelmholtzFluids;
package Propane "Propane"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsPropane},
  helmholtzCoefficients=helmholtzCoefficientsPropane,
  thermalConductivityCoefficients=thermalConductivityCoefficientsPropane,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsPropane,
  surfaceTensionCoefficients=surfaceTensionCoefficientsPropane,
  ancillaryCoefficients=ancillaryCoefficientsPropane,
  fluidLimits=fluidLimitsPropane,
  Density(min=fluidLimitsPropane.DMIN, max=fluidLimitsPropane.DMAX, start=fluidConstantsPropane.molarMass/fluidConstantsPropane.criticalMolarVolume),
  Temperature(min=fluidLimitsPropane.TMIN, max=fluidLimitsPropane.TMAX, start=298.15),
  AbsolutePressure(min=0, max=1000e6, start=101325),
  SpecificEnthalpy(min=fluidLimitsPropane.HMIN, max=fluidLimitsPropane.HMAX, start=(fluidLimitsPropane.HMIN+fluidLimitsPropane.HMAX)/2),
  SpecificEntropy(min=fluidLimitsPropane.SMIN, max=fluidLimitsPropane.SMAX, start=(fluidLimitsPropane.SMIN+fluidLimitsPropane.SMAX)/2));

  final constant FluidConstants
  fluidConstantsPropane(
    chemicalFormula="CH3CH2CH3",
    structureFormula="",
    casRegistryNumber="74-98-6",
    iupacName="",
    molarMass=0.04409562,
    hasCriticalData=true,
       criticalTemperature=369.89,
       criticalPressure=4251200,
       criticalMolarVolume=1/5000,
       HCRIT0=555235.416244571,
       SCRIT0=2051.62550792289,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.084,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.1521,
    triplePointTemperature=85.525,
    triplePointPressure=0.653,
    normalBoilingPoint=231.036,
    meltingPoint=85.525) "Fluid Constants";

  final constant FluidLimits
  fluidLimitsPropane(
    TMIN=fluidConstantsPropane.triplePointTemperature,
    TMAX=650,
    DMIN=Modelica.Constants.small,
    DMAX=909,
    PMIN=Modelica.Constants.small,
    PMAX=1000e6,
    HMIN=-200e3,
    HMAX=+2400e3,
    SMIN=-1400,
    SMAX=9800) "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsPropane(
    idealLog=[
      +3.0,         1.],
    idealPower=[
      -4.970583,        0.;
      +4.293520,        1.],
    idealEinstein=[
      +3.043,          -1.062478;
      +5.874,          -3.344237;
      +9.337,          -5.363757;
      +7.922,         -11.762957],
    residualPoly=[
      +0.42910051E-01,  1.00,  4.,  0.;
       0.17313671E+01,  0.33,  1.,  0.;
      -0.24516524E+01,  0.80,  1.,  0.;
       0.34157466E+00,  0.43,  2.,  0.;
      -0.46047898E+00,  0.90,  2.,  0.],
    residualBwr=[
      -0.66847295E+00,  2.46,  1.,  1.;
       0.20889705E+00,  2.09,  3.,  1.;
       0.19421381E+00,  0.88,  6.,  1.;
      -0.22917851E+00,  1.09,  6.,  1.;
      -0.60405866E+00,  3.25,  2.,  2.;
       0.66680654E-01,  4.62,  3.,  2.],
     residualGauss=[
       0.17534618E-01,  0.76,  1.,  2., 2.,  -0.963,    -2.33,  0.684,  1.283;
       0.33874242E+00,  2.50,  1.,  2., 2.,  -1.977,    -3.47,  0.829,  0.6936;
       0.22228777E+00,  2.75,  1.,  2., 2.,  -1.917,    -3.15,  1.419,  0.788;
      -0.23219062E+00,  3.05,  2.,  2., 2.,  -2.307,    -3.19,  0.817,  0.473;
      -0.92206940E-01,  2.55,  2.,  2., 2.,  -2.546,    -0.92,  1.500,  0.8577;
      -0.47575718E+00,  8.40,  4.,  2., 2.,  -3.28,     -18.8,  1.426,  0.271;
      -0.17486824E-01,  6.75,  1.,  2., 2.,  -14.6,     -547.8, 1.093,  0.948])
  "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsPropane(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=369.85,
    reducingThermalConductivity_0=1,
    lambda_0_num_coeffs=[
      -1.24778E-3,    0;
       8.16371E-3,    1;
       1.99374E-2,    2],
    reducingTemperature_background=369.85,
    reducingMolarVolume_background=1/5000,
    reducingThermalConductivity_background=1,
    lambda_b_coeffs=[
      -3.69500E-2,    0,   1,   0;
       4.82798E-2,    1,   1,   0;
       1.48658E-1,    0,   2,   0;
      -1.35636E-1,    1,   2,   0;
      -1.19986E-1,    0,   3,   0;
       1.17588E-1,    1,   3,   0;
       4.12431E-2,    0,   4,   0;
      -4.36911E-2,    1,   4,   0;
      -4.86905E-3,    0,   5,   0;
       6.16079E-3,    1,   5,   0],
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.875350E-9,
    T_ref=637.68) "Coefficients for the thermal conductivity";

  final constant Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsPropane(
    dynamicViscosityModel=DynamicViscosityModel.VS1,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=0.49748,
    epsilon_kappa=263.88,
    CET=[
      0.141824,   0.50],
    a=[
      0.25104574,   0;
     -0.47271238,   1;
      0.060836515,  3],
    b=[
    -19.572881,       0.00;
     219.73999,      -0.25;
    -1015.3226,      -0.50;
     2471.01251,     -0.75;
    -3375.1717,      -1.00;
     2491.6597,      -1.25;
    -787.26086,      -1.50;
     14.085455,      -2.50;
    -0.34664158,     -5.50],
    reducingTemperature_residual=369.82,
    reducingMolarVolume_residual=1/5000,
    reducingViscosity_residual=1,
    g=[
     0.250053938863E1,      0.0;
     0.215175430074E1,      0.5],
    e=[
     0.359873030195E2,   0.0,    2.00,  0.00,  0;
    -0.180512188564E3,  -1.0,    2.00,  0.00,  0;
     0.877124888223E2,  -2.0,    2.00,  0.00,  0;
    -0.105773052525E3,   0.0,    3.00,  0.00,  0;
     0.205319740877E3,  -1.0,    3.00,  0.00,  0;
    -0.129210932610E3,  -2.0,    3.00,  0.00,  0;
     0.589491587759E2,   0.0,    4.00,  0.00,  0;
    -0.129740033100E3,  -1.0,    4.00,  0.00,  0;
     0.766280419971E2,  -2.0,    4.00,  0.00,  0;
    -0.959407868475E1,   0.0,    5.00,  0.00,  0;
     0.210726986598E2,  -1.0,    5.00,  0.00,  0;
    -0.143971968187E2,  -2.0,    5.00,  0.00,  0;
    -0.161688405374E4,   0.0,    1.00, -1.00,  0],
    nu_po=[
     0.161688405374E4,   0.0,    1,  0,  0],
    de_po=[
     1.,                 0.0,    0,  1,  0;
    -1.,                 0.0,    1,  0,  0])
  "Coefficients for the dynamic viscosity";

  final constant Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsPropane(
    coeffs=[
       0.05666,     1.265;
      -0.005291,    2.265]) "Coefficients for the surface tension";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsPropane(
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -6.7722,       1.0;
       1.6938,       1.5;
      -1.3341,       2.2;
      -3.1876,       4.8;
       0.94937,      6.2],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       1.82205,    0.345;
       0.65802,    0.74;
       0.21109,    2.6;
       0.083973,   7.2],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -2.4887,   0.3785;
      -5.1069,   1.07;
      -12.174,   2.7;
      -30.495,   5.5;
      -52.192,  10.0;
      -134.89,  20.0],
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=85.525,
    p_reducing=0.00017,
    pressureMelting1=[
      -4230000000000.0, 0;
       4230000000001.0, 1.283],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2))
  "Coefficients for the ancillary equations (PS5, DL1, DV3, ML1)";


  annotation (Documentation(info="<html>
These are the coefficients for Propane.

<dl>
<dt> Lemmon, Eric W.; McLinden, M. O. and Wagner, W.</dt>
<dd> <b>Thermodynamic Properties of Propane. III. A Reference Equation of State for Temperatures from the Melting Line to 650 K and Pressures up to 1000 MPa</b><br />
     Journal of Chemical &amp; Engineering Data Vol. 54 Issue 12, Pages 3141-3180 (2009)<br />
     DOI: <a href=\"http://dx.doi.org/10.1021/je900217v\">10.1021/je900217v</a>
</dd>
<dt> Vogel, E.; K&uuml;chenmeister, C.; Bich, E. and Laesecke, A.</dt>
<dd> <b>Reference Correlation of the Viscosity of Propane</b><br />
     Journal of Physical and Chemical Reference Data 27/5, 947-970 (1998)<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.556025\">10.1063/1.556025</a>
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

end Propane;
