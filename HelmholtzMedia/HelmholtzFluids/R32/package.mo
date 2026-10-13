within HelmholtzMedia.HelmholtzFluids;
package R32 "R32 with IIR reference state"
  extends Interfaces.PartialHelmholtzMedium(
    mediumName="R32" "short name",
    fluidConstants={fluidConstantsR32},
    helmholtzCoefficients=helmholtzCoefficientsR32,
    thermalConductivityCoefficients=thermalConductivityCoefficientsR32,
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsR32,
    surfaceTensionCoefficients=surfaceTensionCoefficientsR32,
    ancillaryCoefficients=ancillaryCoefficientsR32,
    fluidLimits=fluidLimitsR32,
    Density(min=fluidLimitsR32.DMIN, max=fluidLimitsR32.DMAX, start=fluidConstantsR32.molarMass/fluidConstantsR32.criticalMolarVolume),
    Temperature(min=fluidLimitsR32.TMIN, max=fluidLimitsR32.TMAX, start=298.15),
    AbsolutePressure(min=0, max=fluidLimitsR32.PMAX, start=101325),
    SpecificEnthalpy(min=fluidLimitsR32.HMIN, max=fluidLimitsR32.HMAX, start=(fluidLimitsR32.HMIN+fluidLimitsR32.HMAX)/2),
    SpecificEntropy(min=fluidLimitsR32.SMIN, max=fluidLimitsR32.SMAX, start=(fluidLimitsR32.SMIN+fluidLimitsR32.SMAX)/2));

  final constant FluidConstants
  fluidConstantsR32(
    chemicalFormula="CH2F2",
    structureFormula="",
    casRegistryNumber="75-10-5",
    iupacName="",
    molarMass=0.052024,
    hasCriticalData=true,
       criticalTemperature=351.255000449,
       criticalPressure=5788581.7687,
       criticalMolarVolume=0.052024/424.435691726 "Denominator=critical density*molecular weight",
       HCRIT0=414146.698495652,
       SCRIT0=1648.64125564321,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=1.978,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.2769,
    triplePointTemperature=136.34,
    triplePointPressure=48 "Search <pressure at triple point>",
    normalBoilingPoint=221.499,
    meltingPoint=221.499 "Search <normal boiling point temperature>") "Fluid Constants";

  final constant EoSLimits
  fluidLimitsR32(
    TMIN=fluidConstantsR32.triplePointTemperature,
    TMAX=435 "upper temperature limit",
    DMIN=Modelica.Constants.eps,
    DMAX=1491,
    PMIN=Modelica.Constants.small,
    PMAX=70e6,
    HMIN=-34e3,
    HMAX=705e3,
    SMIN=-370,
    SMAX=8810) "Helmholtz EoS Limits";

  final constant HelmholtzMedia.Interfaces.Types.HelmholtzCoefficients
  helmholtzCoefficientsR32(
    T_reducing=351.255,
    d_reducing=424.435723,
    idealLog=[
       3.004486E+0,     1.00E0] "CPP <c(i), power of T>-1",
    idealPower=[
      -8.258102267741458,     0.00E0;
       6.353101349786038,     1.00E0] "PH0 <aj, ti for [ai*tau**ti] terms> Different R134a and R32",
    idealEinstein=[
       1.160761, -798.0/351.255;
       2.645151, -4185.0/351.255;
       5.794987, -1806.0/351.255;
       1.129475, -11510.0/351.255] "???",
    residualPoly=[
      0.1046634E+01,   0.250,   1.00;
     -0.5451165E+00,   1.000,   2.00;
     -0.2448595E-02,  -0.250,   5.00;
     -0.4877002E-01,  -1.000,   1.00;
      0.3520158E-01,   2.000,   1.00;
      0.1622750E-02,   2.000,   3.00;
      0.2377225E-04,   0.750,   8.00;
      0.2914900E-01,   0.250,   4.00] "FEQ <a(i),t(i),d(i),l(i)>",
    residualBwr=[
      0.3386203E-02,   18.000,   4.00,   4;
     -0.4202444E-02,   26.000,   4.00,   3;
      0.4782025E-03,   -1.000,   8.00,   1;
     -0.5504323E-02,   25.000,   3.00,   4;
     -0.2418396E-01,    1.750,   5.00,   1;
      0.4209034E+00,    4.000,   1.00,   2;
     -0.4616537E+00,    5.000,   1.00,   2;
     -0.1200513E+01,    1.000,   3.00,   1;
     -0.2591550E+01,    1.500,   1.00,   1;
     -0.1400145E+01,    1.000,   2.00,   1;
      0.8263017E+00,    0.500,   3.00,   1] "See before, why 2 row more?",
   residualGauss=fill(0.0, 0, 9)) "Coefficients of the Helmholtz EoS";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsR32(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=351.255,
    reducingThermalConductivity_0=1.0,
    lambda_0_num_coeffs=[
     0.106548E-1,    0;
    -0.194174E-1,    1;
     0.254295E-1,    2] "<coeff, power in T>",
    reducingTemperature_background=351.255 "<reducing par for T, rho, tcx>",
    reducingMolarVolume_background=1/8150.0846 "See above",
    reducingThermalConductivity_background=1.0 "See above",
    lambda_b_coeffs=[
     0.221878E-01,  0,   1,   0;
    -0.215336E-01,  1,   1,   0;
     0.283523E+00,  0,   2,   0;
    -0.169164E+00,  1,   2,   0;
    -0.297237E+00,  0,   3,   0;
     0.191614E+00,  1,   3,   0;
     0.105727E+00,  0,   4,   0;
    -0.665397E-01,  1,   4,   0;
    -0.123172E-01,  0,   5,   0;
     0.766378E-02,  1,   5,   0],
    xi_0=1.94E-10 "<xi0 (amplitude) [m]>",
    Gamma_0=0.0496 "<gam0 (amplitude) [-]>",
    qd_inverse=5.582925E-10,
    T_ref=526.8825) "Coefficients for the thermal conductivity";

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsR32(
    dynamicViscosityModel=DynamicViscosityModel.VS1,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=0.411 "Lennard-Jones size parameter [nm]",
    epsilon_kappa=290.0 "Lennard-Jones energy parameter [K]",
    reducingTemperature_0=351.255,
    CET=[
      0.577885,    0;
      10.2498,     1;
     -4.95882,     2;
      14.1485,     3;
     -0.816434,    4] "dilute gas, numerator <coeff, power of T/T_red>",
    CET_den=[
      0.896478,    0;
     -0.595706,    1;
      1.0,         2] "dilute gas, denominator <coeff, power of T/T_red>",
    b=[
     -19.572881,    0.0;
      219.73999,   -0.25;
     -1015.3226,   -0.5;
      2471.0125,   -0.75;
     -3375.1717,   -1.0;
      2491.6597,   -1.25;
     -787.26086,   -1.5;
      14.085455,   -2.5;
     -0.34664158,  -5.5] "second viscosity virial coefficient <coeff, power of T/(eps/k)>",
    reducingTemperature_residual=351.255,
    reducingMolarVolume_residual=1/8150.0846,
    e=[
      1.24655,       0.5,   2/3,    0,  0;
      8.85264,       0.5,   5/3,    0,  0;
      0.587282,     -0.5,  14/3,    0,  0;
      2.81507E-06,  -0.5,  44/3,    0,  0;
      4.4106,       -1.5,   8/3,    0,  0] "higher-order term <coeff, power of T/T_red, power of d/d_red, -, ->")
  "Coefficients for the dynamic viscosity";

  final constant HelmholtzMedia.Interfaces.Types.SurfaceTensionCoefficients
  surfaceTensionCoefficientsR32(
    coeffs=[
      0.07147,    1.246]) "Coefficients for the surface tension <sigma0 and n>";

  final constant HelmholtzMedia.Interfaces.Types.AncillaryCoefficients
  ancillaryCoefficientsR32(
    T_crit=351.255,
    p_crit=5782000,
    d_crit=424.435723,
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.74883E+01,   1.0;
       0.19697E+01,   1.5;
      -0.17496E+01,   2.2;
      -0.40224E+01,   4.8;
       0.15209E+01,   6.2] "in <vapor pressure equation>",
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       0.12584E+01,   0.27;
       0.46410E+01,   0.8;
      -0.54870E+01,   1.1;
       0.33115E+01,   1.5;
      -0.61370E+00,   1.8] "in<saturated liquid density equation>",
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -0.22002E+01,   0.336;
      -0.59720E+01,   0.98;
      -0.14571E+02,   2.7;
      -0.42598E+02,   5.7;
       0.42686E+01,   6.5;
      -0.73373E+02,  11.0]) "Coefficients for the ancillary equations in <saturated vapor density equation>";

  annotation (Documentation(info="<html>
These are the coefficients for R32.
The reference state is set to IIR (saturated liquid at 0 &deg;C: h = 200 kJ/kg, s = 1 kJ/(kg K)), the REFPROP default.

<dl>
<dt> Tillner-Roth, R. and Yokozeki, A.</dt>
<dd> <b>An International Standard Equation of State for Difluoromethane (R-32) for Temperatures from the Triple Point at 136.34 K to 435 K and Pressures up to 70 MPa</b><br />
     J. Phys. Chem. Ref. Data, 26:1273-1328, 1997.<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.556002\">10.1063/1.556002</a>
</dd>
<dt> Velliadou, D., Antoniadis, K.D., Assael, M.J. and Huber, M.L.</dt>
<dd> <b>Reference Correlation for the Viscosity of Difluoromethane (R-32) from the Triple Point to 425 K and up to 70 MPa</b><br />
     Int. J. Thermophys., 43:129, 2022.<br />
     DOI: <a href=\"https://doi.org/10.1007/s10765-022-03050-9\">10.1007/s10765-022-03050-9</a>
</dd>
<dt> Perkins, R.A. and Huber, M.L.</dt>
<dd> <b>Thermal conductivity correlation for R-32</b> (2005, unpublished), as implemented in REFPROP 10,
     with the simplified critical enhancement of Olchowy and Sengers (1989).<br />
     Data: DOI <a href=\"https://doi.org/10.1021/je010001m\">10.1021/je010001m</a>
</dd>
<dt> Mulero, A., Cachadi&ntilde;a, I. and Parra, M.I.</dt>
<dd> <b>Recommended Correlations for the Surface Tension of Common Fluids</b><br />
     J. Phys. Chem. Ref. Data, 41:043105, 2012.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.4768782\">10.1063/1.4768782</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end R32;
