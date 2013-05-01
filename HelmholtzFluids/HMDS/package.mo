within HelmholtzMedia.HelmholtzFluids;
package HMDS "HMDS hexamethyldisiloxane (MM)"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsHMDS},
  helmholtzCoefficients=helmholtzCoefficientsHMDS,
  thermalConductivityCoefficients=thermalConductivityCoefficientsHMDS,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsHMDS,
  surfaceTensionCoefficients=surfaceTensionCoefficientsHMDS,
  ancillaryCoefficients=ancillaryCoefficientsHMDS,
  fluidLimits=fluidLimitsHMDS,
  Density(min=fluidLimitsHMDS.DMIN, max=fluidLimitsHMDS.DMAX, start=fluidConstantsHMDS.molarMass/fluidConstantsHMDS.criticalMolarVolume),
  Temperature(min=fluidLimitsHMDS.TMIN, max=fluidLimitsHMDS.TMAX, start=298.15),
  AbsolutePressure(min=0, max=200e6, start=101325),
  SpecificEnthalpy(min=fluidLimitsHMDS.HMIN, max=fluidLimitsHMDS.HMAX, start=(fluidLimitsHMDS.HMIN+fluidLimitsHMDS.HMAX)/2),
  SpecificEntropy(min=fluidLimitsHMDS.SMIN, max=fluidLimitsHMDS.SMAX, start=(fluidLimitsHMDS.SMIN+fluidLimitsHMDS.SMAX)/2));

  final constant Interfaces.PartialHelmholtzMedium.FluidConstants
  fluidConstantsHMDS(
    chemicalFormula="C6H18OSi2",
    structureFormula="",
    casRegistryNumber="107-46-0",
    iupacName="",
    molarMass=0.16237752,
    hasCriticalData=true,
       criticalTemperature=518.75,
       criticalPressure=1939000,
       criticalMolarVolume=1/1589.825,
       HCRIT0=-1,
       SCRIT0=-1,
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

  final constant Interfaces.PartialHelmholtzMedium.FluidLimits
  fluidLimitsHMDS(
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

  final constant Interfaces.PartialHelmholtzMedium.EoS.HelmholtzCoefficients
  helmholtzCoefficientsHMDS(
    idealLog=[
      +3.24680487,         1.],
    idealPower=[
      -5.42495597,        0.;
      4.919495781,         1.],
    idealEinstein=fill(0.0, 0, 2),
    residualPoly=[
       1.01686012,       0.25,    1.0,   0;
      -2.19713029,       1.125,   1.0,   0;
       0.75443188,       1.5,     1.0,   0;
      -0.68003426,       1.375,   2.0,   0;
       0.19082162,       0.25,    3.0,   0;
       0.10530133E-2,    0.875,   7.0,   0],
    residualBwr=[
       0.62845950,       0.625,   2.0,   1;
       0.30903042E-1,    1.75,    5.0,   1;
      -0.83948727,       3.625,   1.0,   2;
      -0.20262381,       3.625,   4.0,   2;
      -0.35131597E-1,   14.5,     3.0,   3;
       0.25902341E-1,   12.0,     4.0,   3],
     residualGauss=fill(0.0, 0, 9)) "Coefficients of the Helmholtz EoS";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsHMDS(
    reducingTemperature_0=1,
    reducingThermalConductivity_0=1,
    lambda_0_num_coeffs=fill(0.0, 0, 2),
    reducingTemperature_residual=1,
    reducingMolarVolume_residual=1,
    reducingThermalConductivity_residual=1,
    lambda_r_coeffs=fill(0.0, 0, 4),
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.875350E-9,
    T_ref=637.68) "Coefficients for the thermal conductivity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsHMDS(
    dynamicViscosityModel=Interfaces.PartialHelmholtzMedium.Types.DynamicViscosityModel.VS1,
    collisionIntegralModel=Interfaces.PartialHelmholtzMedium.Types.CollisionIntegralModel.CI1,
    sigma=1,
    epsilon_kappa=1,
    CET=fill(0.0, 0, 2),
    a=fill(0.0, 0, 2),
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
    reducingTemperature_residual=1,
    reducingMolarVolume_residual=1,
    reducingViscosity_residual=1,
    g=fill(0.0, 0, 2),
    e=fill(0.0, 0, 5),
    nu_po=fill(0.0, 0, 5),
    de_po=fill(0.0, 0, 5)) "Coefficients for the dynamic viscosity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsHMDS(
    coeffs=fill(0.0, 0, 2)) "Coefficients for the surface tension";

  final constant
  Interfaces.PartialHelmholtzMedium.Ancillary.AncillaryCoefficients
  ancillaryCoefficientsHMDS(
    pressureSaturationModel=Interfaces.PartialHelmholtzMedium.Types.PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.86671E+01,   1.0;
       0.11649E+02,   1.5;
      -0.11484E+02,   1.65;
      -0.53256E+01,   4.5],
    densityLiquidModel=Interfaces.PartialHelmholtzMedium.Types.DensityLiquidModel.DL1,
    densityLiquid=[
       0.14533E+02,   0.584;
      -0.49804E+02,   0.80;
       0.83748E+02,   1.02;
      -0.70321E+02,   1.26;
       0.24283E+02,   1.50],
    densityVaporModel=Interfaces.PartialHelmholtzMedium.Types.DensityVaporModel.DV3,
    densityVapor=[
      -0.35719E+01,   0.373;
      -0.14740E+03,   2.15;
       0.40699E+03,   2.6;
      -0.69676E+03,   3.3;
       0.12541E+04,   4.2;
      -0.91199E+03,   4.6])
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for HMDS. 

<dl>
<dt> Bücker, D. and Wagner, W.</dt>
<dd> <b>Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-HMDS and IsoHMDS</b><br>
     Journal of Physical and Chemical Reference Data 35.2, S. 929-1019 (2006)<br>
     DOI: <a href=\"http://dx.doi.org/10.1063/1.1901687\">10.1063/1.1901687</a>
</dd>
<dt> Vogel, Eckhard; Küchenmeister, Cornelia and Bich, Eckard</dt>
<dd> <b>Viscosity correlation for n-HMDS in the fluid region</b><br>
     High Temperatures - High Pressures 31.2, 173-186 (1999)<br>
     DOI: <a href=\"http://dx.doi.org/10.1068/htrt154\">10.1068/htrt154</a>
</dd>
<dt> Perkins, Richard A. et. al.</dt>
<dd> <b>Measurement and Correlation of the Thermal Conductivity of HMDS from 135 K to 600 K at Pressures to 70 MPa</b><br>
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

end HMDS;
