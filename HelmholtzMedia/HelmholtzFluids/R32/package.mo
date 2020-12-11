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
    AbsolutePressure(min=0, max=70e6, start=101325),
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
       criticalTemperature=351.255,
       criticalPressure=5782000,
       criticalMolarVolume=0.052024/424.435723 "Denominator=critical density*molecular weight",
       HCRIT0=0 "How to get this data?",
       SCRIT0=0 "How to get this data?",
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

  final constant FluidLimits
  fluidLimitsR32(
    TMIN=fluidConstantsR32.triplePointTemperature,
    TMAX=435 "upper temperature limit",
    DMIN=Modelica.Constants.small,
    DMAX=1429.3 "=maximum density*molecular weight",
    PMIN=Modelica.Constants.small,
    PMAX=70e6,
    HMIN=0 "How to get this data?",
    HMAX=600e3 "How to get this data?",
    SMIN=0 "How to get this data?",
    SMAX=1e12 "How to get this data?") "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsR32(
    idealLog=[
       3.004486E+0,     1.00E0] "CPP <c(i), power of T>-1",
    idealPower=[
      -3.639189999198519,     0.00E0;
       3.979895667072291,     1.00E0] "PH0 <aj, ti for [ai*tau**ti] terms> Different R134a and R32",
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

  final constant Transport.ThermalConductivityCoefficients
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
    qd_inverse=5.285356E-10,
    T_ref=526.8825) "Coefficients for the thermal conductivity";

    // To be continued
  final constant Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsR32(
    dynamicViscosityModel=DynamicViscosityModel.VS1_alternative,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=0.4098 "<Lennard-Jones coefficient sigma [nm]>",
    epsilon_kappa=289.65 "<Lennard-Jones coefficient epsilon/kappa [K]>",
    CET=[
     0.141824E0, 0.5] "Find VS1, FEQ propane.fld <Chapman-Enskog term>",
    a=[
      0.25104574E+0,  0;
     -0.47271238E+0,  1;
      0.060836515,  3] "<coeff, power of Tstar>",
    b=[
     -19.572881E+0,       0.0;
      219.73999E+0,      -0.25;
     -1015.3226E+0,      -0.5;
      2471.01251E+0,     -0.75;
     -3375.1717E+0,      -1.0;
      2491.6597E+0,      -1.25;
     -787.26086E+0,      -1.5;
      14.085455E+0,      -2.5;
     -0.34664158E+0,     -5.5] "<coeff, power in T* = T/(eps/k)>",
    reducingTemperature_residual=369.82,
    reducingMolarVolume_residual=1/5000,
    reducingViscosity_residual=1 "<reducing parameters for T, rho, eta>",
    g=[
     0,      0.00;
     0,      1.00;
     0,      2.00] "<reducing parameters for T, rho, eta>, next rows",
    e=[
      0.250053938863E+1,   0.0,    0.00,  0.00,  0;
      0.215175430074E+1,   0.5,    0.00,  0.00,  0;
      0.359873030195E+2,   0.0,    2.00,  0.00,  0;
     -0.180512188564E+3,  -1.0,    2.00,  0.00,  0;
      0.877124888223E+2,  -2.0,    2.00,  0.00,  0;
     -0.105773052525E+3,   0.0,    3.00,  0.00,  0;
      0.205319740877E+3,  -1.0,    3.00,  0.00,  0;
     -0.129210932610E+3,  -2.0,    3.00,  0.00,  0;
      0.589491587759E+2,   0.0,    4.00,  0.00,  0;
     -0.129740033100E+3,  -1.0,    4.00,  0.00,  0;
      0.766280419971E+2,  -2.0,    4.00,  0.00,  0;
     -0.959407868475E+1,   0.0,    5.00,  0.00,  0;
      0.210726986598E+2,  -1.0,    5.00,  0.00,  0;
     -0.143971968187E+2,  -2.0,    5.00,  0.00,  0;
     -0.161688405374E+4,   0.0,    1.00, -1.00,  0],
    nu_po=[
     0.161688405374E+4,     0.00,  1.00,  0.00,  0],
    de_po=[
     1.00,                   0.00,  0.00,  1.00,  0;
    -1.00,                   0.00,  1.00,  0.00,  0])
  "Coefficients for the dynamic viscosity";

  final constant Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsR32(
    coeffs=[
      0.07147,    1.246]) "Coefficients for the surface tension <sigma0 and n>";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsR32(
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
These are the coefficients for R134a.
The reference state is set to IIR.

<dl>
<dt> Tillner-Roth, R. and Baehr, H.D.,</dt>
<dd> <b>An international standard formulation of the thermodynamic properties of 1,1,1,2-tetrafluoroethane (HFC-134a) for temperatures from 170 K to 455 K at pressures up to 70 MPa</b><br />
     J. Phys. Chem. Ref. Data, 23:657-729, 1994.<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.555958\">10.1063/1.555958</a>
</dd>
<dt> Huber, M.L., Laesecke, A., and Perkins, R.A.</dt>
<dd> <b>Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a</b><br />
     Ind. Eng. Chem. Res., 42:3163-3178, 2003.<br />
     DOI: <a href=\"http://dx.doi.org/10.1021/ie0300880\">10.1021/ie0300880</a>
</dd>
<dt> Perkins, R.A., Laesecke, A., Howley, J., Ramires, M.L.V., Gurova, A.N., and Cusco, L.</dt>
<dd> <b>Experimental thermal conductivity values for the IUPAC round-robin sample of 1,1,1,2-tetrafluoroethane (R134a)</b><br />
     NISTIR, 2000.
</dd>
<dt> Span, R. and Krauss, R.</dt>
<dd> <b>Properties of R134a</b><br />
     VDI Heat Atlas, Section D2.7, 2010.<br />
     DOI: <a href=\"http://dx.doi.org/10.1007/978-3-540-77877-6_11\">10.1007/978-3-540-77877-6_11</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end R32;
