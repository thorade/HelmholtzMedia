within HelmholtzMedia.HelmholtzFluids;
package Helium "Helium"
extends Interfaces.PartialHelmholtzMedium(
  mediumName="helium" "short name",
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
  SpecificEnthalpy(min=fluidLimitsHelium.HMIN, max=fluidLimitsHelium.HMAX, start=0),
  SpecificEntropy(min=fluidLimitsHelium.SMIN, max=fluidLimitsHelium.SMAX, start=0));

  final constant FluidConstants
  fluidConstantsHelium(
    casRegistryNumber="7440-59-7" "CAS number",
    iupacName="helium-4" "full name",
    structureFormula="He",
    chemicalFormula="He",
    molarMass=0.004002602,
    triplePointTemperature=2.1768,
    normalBoilingPoint=4.2226,
    hasCriticalData=true,
       criticalTemperature=5.1953,
       criticalPressure=227610,
       criticalMolarVolume=1.0/17.3837e3,
       HCRIT0=11229.8913760996,
       SCRIT0=2054.84033693702,
    hasAcentricFactor=true,
       acentricFactor=-0.385,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.0,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    triplePointPressure=5033.5,
    meltingPoint=2.1768) "Fluid Constants";

  final constant FluidLimits
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

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsHelium(
    useLineSearch=true,
    idealLog=[
          1.5,          1],
    idealPower=[
        0.159269591430361000,  0;
        0.476532113807410000,  1],
    idealEinstein=fill(0.0, 0, 2),
    idealCosh=fill(0.0, 0, 2),
    idealSinh=fill(0.0, 0, 2),
    residualPoly=[
          0.014799269,   1.0,     4,  0;
          3.06281562,    0.426,   1,  0;
         -4.25338698,    0.631,   1,  0;
          0.05192797,    0.596,   2,  0;
         -0.165087335,   1.705,   2,  0;
          0.087236897,   0.568,   3,  0],
    residualBwr=[
          2.10653786,    0.9524,  1,  1;
         -0.62835030,    1.471,   1,  2;
         -0.28200301,    1.48,    3,  2;
          1.04234019,    1.393,   2,  1;
         -0.07620555,    3.863,   2,  2;
         -1.35006365,    0.803,   1,  1],
     residualGauss=[
          0.11997252,    3.273,   1,  2, 2,   -8.674,   -8.005,  1.1475,  0.912;
          0.10724500,    0.66,    1,  2, 2,   -4.006,   -1.15,   1.7036,  0.79;
         -0.35374839,    2.629,   1,  2, 2,   -8.1099,  -2.143,  1.6795,  0.90567;
          0.75348862,    1.4379,  2,  2, 2,   -0.1449,  -0.147,  0.9512,  5.1136;
          0.00701871,    3.317,   2,  2, 2,   -0.1784,  -0.154,  4.475,   3.6022;
          0.226283167,   2.3676,  2,  2, 2,   -2.432,   -0.701,  2.7284,  0.6488;
         -0.22464733,    0.7545,  3,  2, 2,   -0.0414,  -0.21,   1.7167,  4.2753;
          0.12413584,    1.353,   2,  2, 2,   -0.421,   -0.134,  1.5237,  2.744;
          0.00901399,    1.982,   2,  2, 2,   -5.8575, -19.256,  0.7649,  0.8736])
  "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsHelium(
    thermalConductivityModel=ThermalConductivityModel.TC0,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK0,
    reducingTemperature_0=10,
    reducingThermalConductivity_0=1,
    lambda_0_num_coeffs=fill(0.0, 0, 2),
    reducingTemperature_background=1,
    reducingMolarVolume_background=1,
    reducingThermalConductivity_background=1,
    lambda_b_coeffs=fill(0.0, 0, 4),
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.875350E-9,
    T_ref=637.68) "Coefficients for the thermal conductivity";

  final constant Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsHelium(
    dynamicViscosityModel=DynamicViscosityModel.VS0,
    collisionIntegralModel=CollisionIntegralModel.CI0,
    sigma=1,
    epsilon_kappa=1,
    CET=fill(0.0, 0, 2),
    a=fill(0.0, 0, 2),
    b=fill(0.0, 0, 2),
    reducingTemperature_residual=1,
    reducingMolarVolume_residual=1,
    reducingViscosity_residual=1,
    g=fill(0.0, 0, 2),
    e=fill(0.0, 0, 5),
    nu_po=fill(0.0, 0, 5),
    de_po=fill(0.0, 0, 5)) "Coefficients for the dynamic viscosity";

  final constant Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsHelium(
    coeffs=[
       0.0004656,   1.04;
       0.001889,    2.468;
      -0.002006,    2.661]) "Coefficients for the surface tension";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsHelium(
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=10,
    p_reducing=1000e3,
    pressureMelting1=[
      -1.7455837,      0.000000;
       1.6979793,      1.555414],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2),
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -3.8357,   1.0;
       1.7062,   1.5;
      -0.71231,  1.25;
       1.0862,   2.8],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       1.0929,   0.286;
       1.6584,   1.2;
      -3.6477,   2.0;
       2.7440,   2.8;
      -2.3859,   6.5],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -1.5789,   0.333;
      -10.749,   1.5;
       17.711,   2.1;
      -15.413,   2.7;
      -14.352,   9.0])
  "Coefficients for the ancillary equations (PS5, DL1, DV3, ML1)";


  annotation (Documentation(info="<html>
These are the coefficients for Helium.

<dl>
<dt> Ortiz-Vega, D.O., Hall, K.R., Arp, V.D., and Lemmon, E.W.</dt>
<dd> <b>Interim equation (final version) for the properties of helium</b><br>
     to be published in J. Phys. Chem. Ref. Data, 2013<br>
     DOI: <a href=\"http://dx.doi.org/\"></a>
</dd>
<dt> McCarty, R.D. and Arp, V.D.</dt>
<dd> <b>A new wide rand equation of state for helium</b><br>
     Adv. Cryo. Eng. 35:1465-1475 (1990)<br>
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
