within HelmholtzFluids;
package Isobutane "Isobutane data, copied from RefProp Isobutane.fld"
  extends HelmholtzFluids.PartialHelmholtzFluid(
  fluidConstants={fluidConstantsIsobutane},
  helmholtzCoefficients=helmholtzCoefficientsIsobutane,
  thermalConductivityCoefficients=thermalConductivityCoefficientsIsobutane,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsIsobutane,
  surfaceTensionCoefficients=surfaceTensionCoefficientsIsobutane,
  ancillaryCoefficients=ancillaryCoefficientsIsobutane,
  fluidLimits=fluidLimitsIsobutane,
  Density(
    min=fluidLimitsIsobutane.DMIN,
    max=fluidLimitsIsobutane.DMAX,
    start=fluidLimitsIsobutane.DNOM),
  Temperature(min=fluidLimitsIsobutane.TMIN, max=fluidLimitsIsobutane.TMAX),
  AbsolutePressure(min=fluidLimitsIsobutane.PMIN, max=fluidLimitsIsobutane.PMAX),
  SpecificEnthalpy(min=fluidLimitsIsobutane.HMIN, max=fluidLimitsIsobutane.HMAX),
  SpecificEntropy(min=fluidLimitsIsobutane.SMIN, max=fluidLimitsIsobutane.SMAX));

  final constant PartialHelmholtzFluid.FluidConstants fluidConstantsIsobutane(
    chemicalFormula="C4H10",
    structureFormula="",
    casRegistryNumber="75-28-5",
    iupacName="",
    molarMass=0.0581222,
    hasCriticalData=true,
       criticalTemperature=407.81,
       criticalPressure=3629000,
       criticalMolarVolume=1/3879.756788,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.132,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.184,
    triplePointTemperature=113.73,
    triplePointPressure=0.0228906605,
    normalBoilingPoint=261.401,
    meltingPoint=113.73) "copied from Isobutane.fld";

  final constant PartialHelmholtzFluid.EosLimits fluidLimitsIsobutane(
    TMIN=fluidConstantsIsobutane.triplePointTemperature,
    TNOM=298.15,
    TMAX=575,
    DMIN=Modelica.Constants.small,
    DNOM=225.5,
    DMAX=800,
    PMIN=Modelica.Constants.small,
    PNOM=101325,
    PMAX=35e6,
    HMIN=-725e3,
    HMAX=+700e3,
    SMIN=-3036,
    SMAX=9283) "Helmholtz EoS Limits";

  final constant PartialHelmholtzFluid.HelmholtzCoefficients helmholtzCoefficientsIsobutane(
    idealLog=[
      +3.05956619,         1.],
    idealPower=[
      +11.60865546,         0.;
      -5.29450411,         1.],
    idealEinstein=[
      +4.94641014,        -0.9512779015;
       4.09475197,        -2.3878958853;
       15.6632824,        -4.3469042691;
       9.73918122,        -10.3688586351],
    residualPoly=[
      +0.20686820727966E+01,    0.5,    1.,   0;
      -0.36400098615204E+01,    1.0,    1.,   0;
       0.51968754427244E+00,    1.5,    1.,   0;
       0.17745845870123E+00,    0.0,    2.,   0;
      -0.12361807851599E+00,    0.5,    3.,   0;
       0.45145314010528E-01,    0.5,    4.,   0;
       0.30476479965980E-01,    0.75,   4.,   0],
    residualBwr=[
      +0.75508387706302E+00,    2.0,    1.,   1;
      -0.85885381015629E+00,    2.5,    1.,   1;
       0.36324009830684E-01,    2.5,    2.,   1;
      -0.19548799450550E-01,    1.5,    7.,   1;
      -0.44452392904960E-02,    1.0,    8.,   1;
       0.46410763666460E-02,    1.5,    8.,   1;
      -0.71444097992825E-01,    4.0,    1.,   2;
      -0.80765060030713E-01,    7.0,    2.,   2;
       0.15560460945053E+00,    3.0,    3.,   2;
       0.20318752160332E-02,    7.0,    3.,   2;
      -0.10624883571689E+00,    3.0,    4.,   2;
       0.39807690546305E-01,    1.0,    5.,   2;
       0.16371431292386E-01,    6.0,    5.,   2;
       0.53212200682628E-03,    0.0,   10.,   2;
      -0.78681561156387E-02,    6.0,    2.,   3;
      -0.30981191888963E-02,   13.0,    6.,   3],
     residualGauss=[
      -0.42276036810382E-01,    2.0,    1.,   2, 2,  -10.,  -150.,  1.16,.85,    0., 0., 0.;
      -0.53001044558079E-02,    0.0,    2.,   2, 2,  -10.,  -200.,  1.13,  1.,   0., 0., 0.])
  "Coefficients of the Helmholtz EoS";

  final constant PartialHelmholtzFluid.ThermalConductivityCoefficients thermalConductivityCoefficientsIsobutane(
    reducingTemperature_0 = 407.85,
    reducingThermalConductivity_0 = 1,
    lambda_0_coeffs=[
     -2.37901E-3,    0;
      1.06601E-2,    1;
      2.15811E-2,    2],
    reducingTemperature_residual=407.85,
    reducingMolarVolume_residual=1/3860,
    reducingThermalConductivity_residual=1,
    lambda_r_coeffs=[
      -4.11789E-2,    0,   1,   0;
       4.76346E-2,    1,   1,   0;
       1.46805E-1,    0,   2,   0;
      -1.28445E-1,    1,   2,   0;
      -1.19190E-1,    0,   3,   0;
       1.07565E-1,    1,   3,   0;
       4.10226E-2,    0,   4,   0;
      -3.85968E-2,    1,   4,   0;
      -4.88704E-3,    0,   5,   0;
       5.20901E-3,    1,   5,   0],
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.657661E-9,
    T_ref=611.73) "Coefficients for the thermal conductivity";

  final constant PartialHelmholtzFluid.DynamicViscosityCoefficients dynamicViscosityCoefficientsIsobutane(
    dynamicViscosityModel=Types.DynamicViscosityModel.VS1,
    sigma=0.46445,
    epsilon_kappa=280.51,
    CET=[
     0.1628213, 0.5],
    a=[
     0.53583008,   0;
    -0.45629630,   1;
     0.049911282,  2],
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
    reducingTemperature_residual=407.817,
    reducingMolarVolume_residual=1/3860,
    reducingViscosity_residual=1,
    g=[
     0.233859774637E1,   0.0;
     0.235255150838E1,   0.5],
    e=[
     0.103511763411E3,   0.0,    2.00,  0.00,  0;
    -0.312670896234E3,  -1.0,    2.00,  0.00,  0;
     0.145253750239E3,  -2.0,    2.00,  0.00,  0;
    -0.210649894193E3,   0.0,    3.00,  0.00,  0;
     0.386269696509E3,  -1.0,    3.00,  0.00,  0;
    -0.214963015527E3,  -2.0,    3.00,  0.00,  0;
     0.112580360920E3,   0.0,    4.00,  0.00,  0;
    -0.223242033154E3,  -1.0,    4.00,  0.00,  0;
     0.119114788598E3,  -2.0,    4.00,  0.00,  0;
    -0.181909745900E2,   0.0,    5.00,  0.00,  0;
     0.360438957232E2,  -1.0,    5.00,  0.00,  0;
    -0.213960184050E2,  -2.0,    5.00,  0.00,  0;
    -0.194037606990E4,   0.0,    1.00, -1.00,  0],
    nu_po=[
     0.194037606990E4,   0.0,    1.00,  0.00,  0],
    de_po=[
      1.,                0.0,    0.00,  1.00,  0;
     -1.,                0.0,    1.00,  0.00,  0])
  "Coefficients for the dynamic viscosity";

    final constant PartialHelmholtzFluid.SurfaceTensionCoefficients surfaceTensionCoefficientsIsobutane(
    coeffs=[
       0.05756,     1.290;
      -0.009554,    2.290]) "Coefficients for the surface tension";

  final constant PartialHelmholtzFluid.AncillaryCoefficients ancillaryCoefficientsIsobutane(
    pressureSaturationModel=Types.PressureSaturationModel.PS5,
    pressureSaturation=[
      -6.85093103,         1.0;
       1.36543198,         1.5;
      -1.32542691,         2.5;
      -2.56190994,         4.5],
    densityLiquidModel=Types.DensityLiquidModel.DL2,
    densityLiquid=[
       2.04025104,          1.065;
       0.850874089,         3.0;
      -0.479052281,         4.0;
       0.348201252,         7.0],
    densityVaporModel=Types.DensityVaporModel.DV6,
    densityVapor=[
      -2.12933323,          1.065;
      -2.93790085,          2.5;
      -0.89441086,          9.5;
      -3.46343707,         13.0])
  "Coefficients for the ancillary equations (PS5, DL2, DV6)";


  annotation (Documentation(info="<html>
These are the coefficients for Isobutane. 
Implementation of the same correlations as in RefProp. 
All data is copied from Isobutane.fld
Units are converted to SI because Modelica uses SI.

<dl>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
<dt> Bücker, D. and Wagner, W.</dt>
<dd> <b>Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Isobutane and Isobutane</b><br>
     Journal of Physical and Chemical Reference Data 35.2, S. 929-1019 (2006)<br>
     DOI: <a href=\"http://dx.doi.org/10.1063/1.1901687\">10.1063/1.1901687</a>
</dd>
<dt> Vogel, Eckhard; Küchenmeister, Cornelia and Bich, Eckard</dt>
<dd> <b>Viscosity correlation for n-Isobutane in the fluid region</b><br>
     High Temperatures - High Pressures 31.2, 173-186 (1999)<br>
     DOI: <a href=\"http://dx.doi.org/10.1068/htrt154\">10.1068/htrt154</a>
</dd>
<dt> Perkins, Richard A. et. al.</dt>
<dd> <b>Measurement and Correlation of the Thermal Conductivity of Isobutane from 135 K to 600 K at Pressures to 70 MPa</b><br>
     Journal of Chemical & Engineering Data 47.5, S. 1263-1271. (2002)<br>
     DOI: <a href=\"http://dx.doi.org/10.1021/je0101202\">10.1021/je0101202</a>
</dd>
</dl>
</html>"));

end Isobutane;
