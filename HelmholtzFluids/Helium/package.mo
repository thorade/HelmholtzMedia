within HelmholtzMedia.HelmholtzFluids;
package Helium "Helium"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsHelium},
  helmholtzCoefficients=helmholtzCoefficientsHelium,
  thermalConductivityCoefficients=thermalConductivityCoefficientsHelium,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsHelium,
  surfaceTensionCoefficients=surfaceTensionCoefficientsHelium,
  ancillaryCoefficients=ancillaryCoefficientsHelium,
  fluidLimits=fluidLimitsHelium,
  Density(min=fluidLimitsHelium.DMIN, max=fluidLimitsHelium.DMAX, start=fluidConstantsHelium.molarMass/fluidConstantsHelium.criticalMolarVolume),
  Temperature(min=fluidLimitsHelium.TMIN, max=fluidLimitsHelium.TMAX, start=298.15),
  AbsolutePressure(min=0, max=1000e6, start=101325),
  SpecificEnthalpy(min=fluidLimitsHelium.HMIN, max=fluidLimitsHelium.HMAX, start=(fluidLimitsHelium.HMIN+fluidLimitsHelium.HMAX)/2),
  SpecificEntropy(min=fluidLimitsHelium.SMIN, max=fluidLimitsHelium.SMAX, start=(fluidLimitsHelium.SMIN+fluidLimitsHelium.SMAX)/2));

  final constant Interfaces.PartialHelmholtzMedium.FluidConstants
  fluidConstantsHelium(
    chemicalFormula="He",
    structureFormula="",
    casRegistryNumber="7440-59-7",
    iupacName="",
    molarMass=0.004002602,
    hasCriticalData=true,
       criticalTemperature=5.1953,
       criticalPressure=227600,
       criticalMolarVolume=1/18.13e3,
       HCRIT0=11338.2625757923,
       SCRIT0=2079.6088606846,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.0,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=-0.385,
    triplePointTemperature=2.1768,
    triplePointPressure=5043.0,
    normalBoilingPoint=4.222,
    meltingPoint=2.1768) "Fluid Constants";

  final constant Interfaces.PartialHelmholtzMedium.FluidLimits
  fluidLimitsHelium(
    TMIN=fluidConstantsHelium.triplePointTemperature,
    TMAX=2000,
    DMIN=Modelica.Constants.small,
    DMAX=565.247,
    PMIN=Modelica.Constants.small,
    PMAX=1000e6,
    HMIN=-6.6877e3,
    HMAX=+12000e3,
    SMIN=-1.8782e3,
    SMAX=120e3) "Helmholtz EoS Limits";

  final constant Interfaces.PartialHelmholtzMedium.EoS.HelmholtzCoefficients
  helmholtzCoefficientsHelium(
    useLineSearch=true,
    idealLog=[
          1.5,          1],
    idealPower=[
        0.187171646656055,  0;
        0.484851136045604,  1],
    idealEinstein=fill(0.0, 0, 2),
    idealCosh=fill(0.0, 0, 2),
    idealSinh=fill(0.0, 0, 2),
    residualPoly=[
        0.9288766E-02,   1.0,    4,  0;
        0.9258069E+00,   0.28,   1,  0;
       -0.1718156E+01,   0.735,  1,  0;
        0.7606137E+00,   0.64,   2,  0;
       -0.1024864E+01,   0.82,   2,  0;
        0.1052455E+00,   1.16,   3,  0],
    residualBwr=[
       -0.1875722E+00,   1.28,   1,  1;
       -0.1287812E+00,   2.0,    1,  2;
       -0.2227619E-02,   0.41,   3,  2;
        0.1823465E+00,   1.33,   2,  1;
       -0.4450014E-01,   4.2,    2,  2;
       -0.8729033E-04,   0.6,    8,  1],
     residualGauss=[
        0.3854320E-01,   3.0,    1,  2, 2,   -1.0833,  -0.0385,  1.9776,  0.6914;
       -0.9585106E+00,   1.0,    1,  2, 2,  -18.3824, -19.8246,  1.6178,  0.8590;
       -0.5454010E-01,   8.2,    1,  2, 2,   -5.0573,  -9.3799,  0.4371,  0.8787;
       -0.3687260E-01,   1.0,    2,  2, 2,   -0.2832,  -0.8073,  0.5355,  2.7182;
       -0.1021851E-02,   2.71,   2,  2, 2,   -6.0582,  -0.0310,  0.7777,  2.0301;
        0.6166348E-01,   1.0,    2,  2, 2,   -0.2444,  -0.0061,  0.4832,  0.8900;
        0.2493437E-01,   1.0,    3,  2, 2,   -0.0539,  -0.3581,  0.8162,  1.1790;
       -0.8127424E-02,   2.0,    3,  2, 2,   -0.1850,  -0.7518,  1.2896,  0.5680;
       -0.8233032E-02,   1.0,    2,  2, 2,   -0.5941,  -7.4629,  0.3577,  1.6412])
  "Coefficients of the Helmholtz EoS";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsHelium(
    reducingTemperature_0=425.16,
    reducingThermalConductivity_0=1,
    lambda_0_num_coeffs=[
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

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsHelium(
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

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsHelium(
    coeffs=[
      0.05418,    1.26]) "Coefficients for the surface tension";

  final constant
  Interfaces.PartialHelmholtzMedium.Ancillary.AncillaryCoefficients
  ancillaryCoefficientsHelium(
    pressureSaturationModel=Interfaces.PartialHelmholtzMedium.Types.PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.399865E+01, 1.0;
       0.870145E+00, 1.5;
       0.171451E+00, 1.85;
       0.120927E+01, 2.7],
    densityLiquidModel=Interfaces.PartialHelmholtzMedium.Types.DensityLiquidModel.DL2,
    densityLiquid=[
       0.140808E+01, 1.17;
      -0.543843E+00, 7.0;
       0.177220E+01, 15.0;
      -0.344056E+01, 20.0],
    densityVaporModel=Interfaces.PartialHelmholtzMedium.Types.DensityVaporModel.DV3,
    densityVapor=[
      -0.126074E+01, 0.263;
      -0.363425E+01, 1.04;
      -0.487998E+01, 3.25;
      -0.130581E+02, 8.5],
    pressureMeltingModel=Interfaces.PartialHelmholtzMedium.Types.PressureMeltingModel.ML1,
    T_reducing=1,
    p_reducing=1000e3,
    pressureMelting1=[
      -1.7455837, 0.000000;
       1.6979793, 1.555414],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2))
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for Helium. 
Warning: The transport properties are those of n-Butane!

<dl>
<dt> Ortiz-Vega, D.O., Hall, K.R., Arp, V.D., and Lemmon, E.W.</dt>
<dd> <b>Interim equation for the properties of helium</b><br>
     to be published in Int. J. Thermophys.<br>
     DOI: <a href=\"http://dx.doi.org/\"></a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Helium;
