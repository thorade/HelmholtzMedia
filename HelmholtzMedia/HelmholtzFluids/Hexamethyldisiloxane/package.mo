within HelmholtzMedia.HelmholtzFluids;
package Hexamethyldisiloxane "hexamethyldisiloxane (MM)"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsMM},
  helmholtzCoefficients=helmholtzCoefficientsMM,
  thermalConductivityCoefficients=thermalConductivityCoefficientsMM,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsMM,
  surfaceTensionCoefficients=surfaceTensionCoefficientsMM,
  ancillaryCoefficients=ancillaryCoefficientsMM,
  fluidLimits=fluidLimitsMM,
  Density(min=fluidLimitsMM.DMIN, max=fluidLimitsMM.DMAX, start=fluidConstantsMM.molarMass/fluidConstantsMM.criticalMolarVolume),
  Temperature(min=fluidLimitsMM.TMIN, max=fluidLimitsMM.TMAX, start=298.15),
  AbsolutePressure(min=0, max=200e6, start=101325),
  SpecificEnthalpy(min=fluidLimitsMM.HMIN, max=fluidLimitsMM.HMAX, start=0),
  SpecificEntropy(min=fluidLimitsMM.SMIN, max=fluidLimitsMM.SMAX, start=0));

  final constant FluidConstants
  fluidConstantsMM(
    chemicalFormula="C6H18OSi2",
    structureFormula="",
    casRegistryNumber="107-46-0",
    iupacName="",
    molarMass=0.16237752,
    hasCriticalData=true,
       criticalTemperature=518.69997204,
       criticalPressure=1939390,
       criticalMolarVolume=1/1874.67076,
       HCRIT0=364989.364617502,
       SCRIT0=806.290249664929,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=false,
       dipoleMoment=-1,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.418,
    triplePointTemperature=204.93,
    triplePointPressure=0002.69,
    normalBoilingPoint=373.401,
    meltingPoint=204.93) "Fluid Constants";

  final constant FluidLimits
  fluidLimitsMM(
    TMIN=273,
    TMAX=673,
    DMIN=Modelica.Constants.small,
    DMAX=844,
    PMIN=Modelica.Constants.small,
    PMAX=30e6,
    HMIN=-100e3,
    HMAX=+1300e3,
    SMIN=-465,
    SMAX=8100) "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsMM(
    idealLog=[
      +51.894/8.314472-1,           1],
    idealPower=[
      +27.190399518670160,                       0;
      -7.3982398876225090,                       1;
      -741.34e-3/8.314472 /2 *518.69997204^1,   -1;
      +416.10e-6/8.314472 /6 *518.69997204^2,   -2;
      -70.000e-9/8.314472 /12*518.69997204^3,   -3],
    idealEinstein=fill(0.0, 0, 2),
    residualPoly=[
       1.01686012,       0.25,    1,   0;
      -2.19713029,       1.125,   1,   0;
       0.75443188,       1.5,     1,   0;
      -0.68003426,       1.375,   2,   0;
       0.19082162,       0.25,    3,   0;
       0.10530133E-2,    0.875,   7,   0],
    residualBwr=[
       0.62845950,       0.625,   2,   1;
       0.30903042E-1,    1.75,    5,   1;
      -0.83948727,       3.625,   1,   2;
      -0.20262381,       3.625,   4,   2;
      -0.35131597E-1,   14.5,     3,   3;
       0.25902341E-1,   12.0,     4,   3],
     residualGauss=fill(0.0, 0, 9)) "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsMM(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=1,
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
  dynamicViscosityCoefficientsMM(
    dynamicViscosityModel=DynamicViscosityModel.VS1,
    collisionIntegralModel=CollisionIntegralModel.CI1,
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
  surfaceTensionCoefficientsMM(
    coeffs=fill(0.0, 0, 2)) "Coefficients for the surface tension";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsMM(
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.86671E+01,   1.0;
       0.11649E+02,   1.5;
      -0.11484E+02,   1.65;
      -0.53256E+01,   4.5],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       0.14533E+02,   0.584;
      -0.49804E+02,   0.80;
       0.83748E+02,   1.02;
      -0.70321E+02,   1.26;
       0.24283E+02,   1.50],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -0.35719E+01,   0.373;
      -0.14740E+03,   2.15;
       0.40699E+03,   2.6;
      -0.69676E+03,   3.3;
       0.12541E+04,   4.2;
      -0.91199E+03,   4.6])
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for MM.

<dl>
<dt> Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W.</dt>
<dd> <b>Multiparameter Equations of State for Selected Siloxanes</b><br>
     Fluid Phase Equilibria, 244:193-211, (2006)<br>
     DOI: <a href=\"http://dx.doi.org/10.1016/j.fluid.2006.04.015\">10.1016/j.fluid.2006.04.015</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Hexamethyldisiloxane;
