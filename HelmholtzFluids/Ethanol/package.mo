within HelmholtzMedia.HelmholtzFluids;
package Ethanol "Ethanol"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsEthanol},
  helmholtzCoefficients=helmholtzCoefficientsEthanol,
  thermalConductivityCoefficients=thermalConductivityCoefficientsEthanol,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsEthanol,
  surfaceTensionCoefficients=surfaceTensionCoefficientsEthanol,
  ancillaryCoefficients=ancillaryCoefficientsEthanol,
  fluidLimits=fluidLimitsEthanol,
  Density(min=fluidLimitsEthanol.DMIN, max=fluidLimitsEthanol.DMAX),
  Temperature(min=fluidLimitsEthanol.TMIN, max=fluidLimitsEthanol.TMAX),
  AbsolutePressure(min=0, max=200e6),
  SpecificEnthalpy(min=fluidLimitsEthanol.HMIN, max=fluidLimitsEthanol.HMAX),
  SpecificEntropy(min=fluidLimitsEthanol.SMIN, max=fluidLimitsEthanol.SMAX));

  final constant Interfaces.PartialHelmholtzMedium.FluidConstants
  fluidConstantsEthanol(
    chemicalFormula="C2H6O",
    structureFormula="",
    casRegistryNumber="64-17-5",
    iupacName="",
    molarMass=0.04606844,
    hasCriticalData=true,
       criticalTemperature=513.9,
       criticalPressure=6148000,
       criticalMolarVolume=0.04606844/276,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=1.6909,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.644,
    triplePointTemperature=159,
    triplePointPressure=0.00043,
    normalBoilingPoint=351.39,
    meltingPoint=159) "Fluid Constants";

  final constant Interfaces.PartialHelmholtzMedium.FluidLimits
  fluidLimitsEthanol(
    TMIN=fluidConstantsEthanol.triplePointTemperature,
    TMAX=650,
    DMIN=Modelica.Constants.small,
    DMAX=893.73,
    PMIN=Modelica.Constants.small,
    PMAX=280e6,
    HMIN=-100e3,
    HMAX=+1300e3,
    SMIN=-465,
    SMAX=8100) "Helmholtz EoS Limits";

  final constant Interfaces.PartialHelmholtzMedium.HelmholtzCoefficients
  helmholtzCoefficientsEthanol(
    idealLog=[
      +5.41129,         1],
    idealPower=[
      -1,        0;
      +1,        1],
    idealEinstein=[
      +1.95989, -694/513.9;
      +7.60084, -1549/513.9;
      +3.89583, -2911/513.9;
      +4.23238, -4659/513.9],
    residualPoly=[
      +0.114008942201E+2,  -0.5,   1,   0;
      -0.395227128302E+2,   0.0,   1,   0;
       0.413063408370E+2,   0.5,   1,   0;
      -0.188892923721E+2,   1.5,   1,   0;
       0.472310314140E+1,   2.0,   1,   0;
      -0.778322827052E-2,   5.0,   1,   0;
       0.171707850032E+0,  -0.5,   2,   0;
      -0.153758307602E+1,   1.0,   2,   0;
       0.142405508571E+1,   2.0,   2,   0;
       0.132732097050E+0,   0.0,   3,   0;
      -0.114231649761E+0,   2.5,   3,   0;
       0.327686088736E-5,   6.0,   6,   0;
       0.495699527725E-3,   2.0,   7,   0;
      -0.701090149558E-4,   2.0,   8,   0;
      -0.225019381648E-5,   4.0,   8,   0],
    residualBwr=[
      -0.255406026981E+0,   5.0,   1,   2;
      -0.632036870646E-1,   3.0,   3,   2;
      -0.314882729522E-1,   7.0,   3,   2;
       0.256187828185E-1,   5.5,   6,   2;
      -0.308694499382E-1,   4.0,   7,   2;
       0.722046283076E-2,   1.0,   8,   2;
       0.299286406225E-2,  22.0,   2,   4;
       0.972795913095E-3,  23.0,   7,   4],
     residualGauss=fill(0.0, 0, 12)) "Coefficients of the Helmholtz EoS";

  final constant
  Interfaces.PartialHelmholtzMedium.ThermalConductivityCoefficients
  thermalConductivityCoefficientsEthanol(
    reducingTemperature_0=425.16,
    reducingThermalConductivity_0=1,
    lambda_0_coeffs=[
     1.62676E-03,    0;
     9.75703E-04,    1;
     2.89887E-02,    2],
    reducingTemperature_residual=425.16,
    reducingMolarVolume_residual=1/3920,
    reducingThermalConductivity_residual=1,
    lambda_r_coeffs=[
    -3.04337E-2,    0,   1,   0;
     4.18357E-2,    1,   1,   0;
     1.65820E-1,    0,   2,   0;
    -1.47163E-1,    1,   2,   0;
    -1.48144E-1,    0,   3,   0;
     1.33542E-1,    1,   3,   0;
     5.25500E-2,    0,   4,   0;
    -4.85489E-2,    1,   4,   0;
    -6.29367E-3,    0,   5,   0;
     6.44307E-3,    1,   5,   0],
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.875350E-9,
    T_ref=637.68) "Coefficients for the thermal conductivity";

  final constant Interfaces.PartialHelmholtzMedium.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsEthanol(
    dynamicViscosityModel=Interfaces.PartialHelmholtzMedium.Types.DynamicViscosityModel.VS1,
    collisionIntegralModel=Interfaces.PartialHelmholtzMedium.Types.CollisionIntegralModel.CI1,
    sigma=0.57335,
    epsilon_kappa=280.51,
    CET=[
     0.1628213, 0.5],
    a=[
     0.17067154,    0;
    -0.48879666,    1;
     0.039038856,   2],
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
    reducingTemperature_residual=425.125,
    reducingMolarVolume_residual=1/3920,
    reducingViscosity_residual=1,
    g=[
     2.30873963359,      0.0;
     2.03404037254,      0.5],
    e=[
    -54.7737770846,      0.0,    2,  0,  0;
     58.0898623034,     -1.0,    2,  0,  0;
     0,                 -2.0,    2,  0,  0;
     35.2658446259,      0.0,    3,  0,  0;
    -39.6682203832,     -1.0,    3,  0,  0;
     0,                 -2.0,    3,  0,  0;
    -1.83729542151,      0.0,    4,  0,  0;
     0,                 -1.0,    4,  0,  0;
     0,                 -2.0,    4,  0,  0;
    -0.833262985358,     0.0,    5,  0,  0;
     1.93837020663,     -1.0,    5,  0,  0;
     0,                 -2.0,    5,  0,  0;
    -188.075903903,      0.0,    1, -1,  0],
    nu_po=[
     188.075903903,      0.0,    1,  0,  0],
    de_po=[
     1.,                 0.0,    0,  1,  0;
    -1.,                 0.0,    1,  0,  0])
  "Coefficients for the dynamic viscosity";

  final constant Interfaces.PartialHelmholtzMedium.SurfaceTensionCoefficients
  surfaceTensionCoefficientsEthanol(
    coeffs=[
      0.065,    1.26]) "Coefficients for the surface tension";

  final constant Interfaces.PartialHelmholtzMedium.AncillaryCoefficients
  ancillaryCoefficientsEthanol(
    pressureSaturationModel=Interfaces.PartialHelmholtzMedium.Types.PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.81829E+01,   1.0;
      -0.62767E+00,   1.5;
      -0.33289E+01,   3.0;
      -0.35278E+01,   5.6;
       0.93103E+01,   7.0],
    densityLiquidModel=Interfaces.PartialHelmholtzMedium.Types.DensityLiquidModel.DL1,
    densityLiquid=[
       0.11818E+01,   0.098;
      -0.36120E+01,   0.22;
       0.54325E+01,   0.35;
      -0.47789E+00,   0.7;
      -0.17766E-01,   2.0],
    densityVaporModel=Interfaces.PartialHelmholtzMedium.Types.DensityVaporModel.DV3,
    densityVapor=[
      -0.93315E+00,   0.09;
      -0.40761E+02,   1.07;
       0.63250E+02,   1.3;
      -0.45195E+02,   1.7;
       0.15114E+01,   4.0;
      -0.56666E+02,   5.0])
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for Ethanol. 

<dl>
<dt> Bücker, D. and Wagner, W.</dt>
<dd> <b>Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Ethanol and IsoEthanol</b><br>
     Journal of Physical and Chemical Reference Data 35.2, S. 929-1019 (2006)<br>
     DOI: <a href=\"http://dx.doi.org/10.1063/1.1901687\">10.1063/1.1901687</a>
</dd>
<dt> Vogel, Eckhard; Küchenmeister, Cornelia and Bich, Eckard</dt>
<dd> <b>Viscosity correlation for n-Ethanol in the fluid region</b><br>
     High Temperatures - High Pressures 31.2, 173-186 (1999)<br>
     DOI: <a href=\"http://dx.doi.org/10.1068/htrt154\">10.1068/htrt154</a>
</dd>
<dt> Perkins, Richard A. et. al.</dt>
<dd> <b>Measurement and Correlation of the Thermal Conductivity of Ethanol from 135 K to 600 K at Pressures to 70 MPa</b><br>
     Journal of Chemical & Engineering Data 47.5, S. 1263-1271. (2002)<br>
     DOI: <a href=\"http://dx.doi.org/10.1021/je0101202\">10.1021/je0101202</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Ethanol;
