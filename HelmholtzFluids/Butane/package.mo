within HelmholtzMedia.HelmholtzFluids;
package Butane "Butane data, copied from RefProp Butane.fld"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsButane},
  helmholtzCoefficients=helmholtzCoefficientsButane,
  thermalConductivityCoefficients=thermalConductivityCoefficientsButane,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsButane,
  surfaceTensionCoefficients=surfaceTensionCoefficientsButane,
  ancillaryCoefficients=ancillaryCoefficientsButane,
  fluidLimits=fluidLimitsButane,
  Density(min=fluidLimitsButane.DMIN, max=fluidLimitsButane.DMAX),
  Temperature(min=fluidLimitsButane.TMIN, max=fluidLimitsButane.TMAX),
  AbsolutePressure(min=0, max=200e6),
  SpecificEnthalpy(min=fluidLimitsButane.HMIN, max=fluidLimitsButane.HMAX),
  SpecificEntropy(min=fluidLimitsButane.SMIN, max=fluidLimitsButane.SMAX));

  final constant Interfaces.PartialHelmholtzMedium.FluidConstants
  fluidConstantsButane(
    chemicalFormula="C4H10",
    structureFormula="",
    casRegistryNumber="106-97-8",
    iupacName="",
    molarMass=0.0581222,
    hasCriticalData=true,
       criticalTemperature=425.125,
       criticalPressure=3796000,
       criticalMolarVolume=0.0581222/228,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.05,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.201,
    triplePointTemperature=134.895,
    triplePointPressure=0.653,
    normalBoilingPoint=272.660,
    meltingPoint=134.912) "copied from Butane.fld";

  final constant Interfaces.PartialHelmholtzMedium.FluidLimits
  fluidLimitsButane(
    TMIN=fluidConstantsButane.triplePointTemperature,
    TMAX=575,
    DMIN=Modelica.Constants.small,
    DMAX=800,
    PMIN=Modelica.Constants.small,
    PMAX=200e6,
    HMIN=-725e3,
    HMAX=+700e3,
    SMIN=-3036,
    SMAX=9283) "Helmholtz EoS Limits";

  final constant Interfaces.PartialHelmholtzMedium.HelmholtzCoefficients
  helmholtzCoefficientsButane(
    idealLog=[
      +3.24680487,         1.],
    idealPower=[
      +12.54882924,        0.;
      -5.46976878,         1.],
    idealEinstein=[
      +5.54913289,        -0.7748404445;
      +11.4648996,        -3.3406025522;
      +7.59987584,        -4.9705130961;
      +9.66033239,        -9.9755537783],
    residualPoly=[
      +0.25536998241635E+01,    0.5,    1.,   0;
      -0.44585951806696E+01,    1.0,    1.,   0;
      +0.82425886369063E+00,    1.5,    1.,   0;
      +0.11215007011442E+00,    0.0,    2.,   0;
      -0.35910933680333E-01,    0.5,    3.,   0;
      +0.16790508518103E-01,    0.5,    4.,   0;
      +0.32734072508724E-01,    0.75,   4.,   0],
    residualBwr=[
      +0.95571232982005E+00,    2.0,    1.,   1;
      -0.10003385753419E+01,    2.5,    1.,   1;
      +0.85581548803855E-01,    2.5,    2.,   1;
      -0.25147918369616E-01,    1.5,    7.,   1;
      -0.15202958578918E-02,    1.0,    8.,   1;
      +0.47060682326420E-02,    1.5,    8.,   1;
      -0.97845414174006E-01,    4.0,    1.,   2;
      -0.48317904158760E-01,    7.0,    2.,   2;
      +0.17841271865468E+00,    3.0,    3.,   2;
      +0.18173836739334E-01,    7.0,    3.,   2;
      -0.11399068074953E+00,    3.0,    4.,   2;
      +0.19329896666669E-01,    1.0,    5.,   2;
      +0.11575877401010E-02,    6.0,    5.,   2;
      +0.15253808698116E-03,    0.0,   10.,   2;
      -0.43688558458471E-01,    6.0,    2.,   3;
      -0.82403190629989E-02,   13.0,    6.,   3],
     residualGauss=[
       -0.28390056949441E-01,   2.0,    1., 2, 2,  -10.,  -150.,  1.16,.85,       0., 0., 0.;
       +0.14904666224681E-02,   0.0,    2., 2, 2,  -10.,  -200.,  1.13,  1.,      0., 0., 0.])
  "Coefficients of the Helmholtz EoS";

  final constant
  Interfaces.PartialHelmholtzMedium.ThermalConductivityCoefficients
  thermalConductivityCoefficientsButane(
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
  dynamicViscosityCoefficientsButane(
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
  surfaceTensionCoefficientsButane(
    coeffs=[
      0.05418,    1.26]) "Coefficients for the surface tension";

  final constant Interfaces.PartialHelmholtzMedium.AncillaryCoefficients
  ancillaryCoefficientsButane(
    pressureSaturationModel=Interfaces.PartialHelmholtzMedium.Types.PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.71897E+01,   1.0;
       0.26122E+01,   1.5;
      -0.21729E+01,   2.0;
      -0.27230E+01,   4.5],
    densityLiquidModel=Interfaces.PartialHelmholtzMedium.Types.DensityLiquidModel.DL1,
    densityLiquid=[
       0.52341E+01,   0.44;
      -0.62011E+01,   0.60;
       0.36063E+01,   0.76;
       0.22137E+00,   5.00],
    densityVaporModel=Interfaces.PartialHelmholtzMedium.Types.DensityVaporModel.DV3,
    densityVapor=[
      -0.27390E+01,   0.391;
      -0.57347E+01,   1.14;
      -0.16408E+02,   3.0;
      -0.46986E+02,   6.5;
      -0.10090E+03,  14.0])
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for Butane. 
Implementation of the same correlations as in RefProp. 
All data is copied from butane.fld
Units are converted to SI because Modelica uses SI.

<dl>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
<dt> Bücker, D. and Wagner, W.</dt>
<dd> <b>Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Butane and Isobutane</b><br>
     Journal of Physical and Chemical Reference Data 35.2, S. 929-1019 (2006)<br>
     DOI: <a href=\"http://dx.doi.org/10.1063/1.1901687\">10.1063/1.1901687</a>
</dd>
<dt> Vogel, Eckhard; Küchenmeister, Cornelia and Bich, Eckard</dt>
<dd> <b>Viscosity correlation for n-butane in the fluid region</b><br>
     High Temperatures - High Pressures 31.2, 173-186 (1999)<br>
     DOI: <a href=\"http://dx.doi.org/10.1068/htrt154\">10.1068/htrt154</a>
</dd>
<dt> Perkins, Richard A. et. al.</dt>
<dd> <b>Measurement and Correlation of the Thermal Conductivity of Butane from 135 K to 600 K at Pressures to 70 MPa</b><br>
     Journal of Chemical & Engineering Data 47.5, S. 1263-1271. (2002)<br>
     DOI: <a href=\"http://dx.doi.org/10.1021/je0101202\">10.1021/je0101202</a>
</dd>
</dl>
</html>"));

end Butane;
