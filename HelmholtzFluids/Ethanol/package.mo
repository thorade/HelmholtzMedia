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
  Temperature(min=fluidLimitsEthanol.TMIN, max=fluidLimitsEthanol.TMAX, start=298.15),
  AbsolutePressure(min=0, max=280e6, start=101325),
  SpecificEnthalpy(min=fluidLimitsEthanol.HMIN, max=fluidLimitsEthanol.HMAX, start=(fluidLimitsEthanol.HMIN+fluidLimitsEthanol.HMAX)/2),
  SpecificEntropy(min=fluidLimitsEthanol.SMIN, max=fluidLimitsEthanol.SMAX, start=(fluidLimitsEthanol.SMIN+fluidLimitsEthanol.SMAX)/2));

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

  final constant Interfaces.PartialHelmholtzMedium.EoS.HelmholtzCoefficients
  helmholtzCoefficientsEthanol(
    idealLog=[
      +5.41129,         1],
    idealPower=[
      -12.899375,        0;
      +8.6222292,        1],
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
  Interfaces.PartialHelmholtzMedium.Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsEthanol(
    reducingTemperature_0=513.9,
    reducingThermalConductivity_0=1,
    lambda_0_coeffs=[
     0.123120E-01,    0;
    -0.153612E-01,    1;
     0.426611E-01,    2],
    reducingTemperature_residual=513.9,
    reducingMolarVolume_residual=1/5991,
    reducingThermalConductivity_residual=1,
    lambda_r_coeffs=[
     0.266894E-1,        0.0,   1.0,   0.0;
     0.0,                1.0,   1.0,   0.0;
    -0.482953E-01,       0.0,   2.0,   0.0;
     0.414022E-01,       1.0,   2.0,   0.0;
     0.172939E-01,       0.0,   3.0,   0.0;
    -0.977825E-02,       1.0,   3.0,   0.0;
     0.0,                0.0,   4.0,   0.0;
     0.0,                1.0,   4.0,   0.0;
     0.0,                0.0,   5.0,   0.0;
     0.0,                1.0,   5.0,   0.0],
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.875350E-9,
    T_ref=637.68) "Coefficients for the thermal conductivity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsEthanol(
    dynamicViscosityModel=Interfaces.PartialHelmholtzMedium.Types.DynamicViscosityModel.VS1,
    collisionIntegralModel=Interfaces.PartialHelmholtzMedium.Types.CollisionIntegralModel.CI1,
    sigma=0.453,
    epsilon_kappa=362.6,
    CET=[
     0, 0.5;
    -1.03116E+0,    0;
     3.48379E-2,    1;
    -6.50264E-6,   2],
    a=[
     0,    0;
     0,    1],
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
    reducingTemperature_residual=513.9,
    reducingMolarVolume_residual=1/5991,
    reducingViscosity_residual=1E3,
    g=[
     -3.38264465E+00,      0.0;
      1.27568864E+01,      0.5],
    e=[
     1.31194057E-01,   0.0,   2.00,  0.00,  0;
    -8.05700894E-02,   0.0,   3.00,  0.00,  0;
    -3.82240694E-01,  -1.00,  2.00,  0.00,  0;
     1.53811778E-01,  -1.00,  3.00,  0.00,  0;
     0.0,             -2.00,  2.00,  0.00,  0;
    -1.10578307E-01,  -2.00,  3.00,  0.00,  0;
    -2.37222995E+01,   0.00,  1.00, -1.00,  0],
    nu_po=[
     2.37222995E+01,   0.00,  1.00,  0.00,  0],
    de_po=[
     1.,                 0.0,    0,  1,  0;
    -1.,                 0.0,    1,  0,  0])
  "Coefficients for the dynamic viscosity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsEthanol(
    coeffs=[
      0.065,    1.26]) "Coefficients for the surface tension";

  final constant
  Interfaces.PartialHelmholtzMedium.Ancillary.AncillaryCoefficients
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
<dt> Dillon, H. E. & Penoncello, S. G.</dt>
<dd> <b>A Fundamental Equation for Calculation of the Thermodynamic Properties of Ethanol</b><br>
     International Journal of Thermophysics 25, 321-335 (2004)<br>
     DOI: <a href=\"http://dx.doi.org/10.1023/B:IJOT.0000028470.49774.14\">10.1023/B:IJOT.0000028470.49774.14</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Ethanol;
