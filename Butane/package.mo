within HelmholtzFluids;
package Butane "Butane data, copied from RefProp Butane.fld"
  extends HelmholtzFluids.PartialHelmholtzFluid(
  fluidConstants={fluidConstantsButane},
  helmholtzCoefficients=helmholtzCoefficientsButane,
  ancillaryCoefficients=ancillaryCoefficientsButane,
  thermalConductivityCoefficients=thermalConductivityCoefficientsButane,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsButane,
  surfaceTensionCoefficients=surfaceTensionCoefficientsButane,
  fluidLimits=fluidLimitsButane,
  Density(
    min=fluidLimitsButane.DMIN,
    max=fluidLimitsButane.DMAX,
    start=fluidLimitsButane.DNOM),
  Temperature(min=fluidLimitsButane.TMIN, max=fluidLimitsButane.TMAX),
  AbsolutePressure(min=fluidLimitsButane.PMIN, max=fluidLimitsButane.PMAX),
  SpecificEnthalpy(min=fluidLimitsButane.HMIN, max=fluidLimitsButane.HMAX),
  SpecificEntropy(min=fluidLimitsButane.SMIN, max=fluidLimitsButane.SMAX));

  final constant PartialHelmholtzFluid.FluidConstants fluidConstantsButane(
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

  final constant PartialHelmholtzFluid.EosLimits fluidLimitsButane(
    TMIN=fluidConstantsButane.triplePointTemperature,
    TNOM=298.15,
    TMAX=575,
    DMIN=Modelica.Constants.small,
    DNOM=228,
    DMAX=800,
    PMIN=Modelica.Constants.small,
    PNOM=101325,
    PMAX=200e6,
    HMIN=-725e3,
    HMAX=+700e3,
    SMIN=-3036,
    SMAX=9283) "Helmholtz EoS Limits";

  final constant PartialHelmholtzFluid.HelmholtzCoefficients helmholtzCoefficientsButane(
    n_ideal={12.54882924,-5.46976878,3.24680487,5.54913289,11.4648996,7.59987584,
        9.66033239},
    Theta={0,0,0,0.7748404445,3.3406025522,4.9705130961,9.9755537783},
    n_residual={0.25536998241635E+01,-0.4458595180669E+01,0.82425886369063,0.11215007011442,
        -0.3591093368033E-01,0.1679050851810E-01,0.3273407250872E-01,0.95571232982005,
        -0.1000338575341E+01,0.8558154880385E-01,-0.2514791836961E-01,-0.15202958578918E-02,
        0.47060682326420E-02,-0.97845414174006E-01,-0.48317904158760E-01,0.17841271865468,
        0.18173836739334E-01,-0.11399068074953,0.19329896666669E-01,0.11575877401010E-02,
        0.15253808698116E-03,-0.43688558458471E-01,-0.82403190629989E-02,-0.28390056949441E-01,
        0.14904666224681E-02},
    c={0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,0,0},
    d={1,1,1,2,3,4,4,1,1,2,7,8,8,1,2,3,3,4,5,5,10,2,6,1,2},
    t={0.50,1.00,1.50,0.00,0.50,0.50,0.75,2.00,2.50,2.50,1.50,1.00,1.50,4.00,7.00,
        3.00,7.00,3.00,1.00,6.00,0.00,6.00,13.00,2.00,0.00},
    crit_eta={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,10},
    crit_beta={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,150,200},
    crit_epsilon={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.85,1.00},
    crit_gamma={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.16,1.13})
  "Coefficients of the Helmholtz EoS";

  final constant PartialHelmholtzFluid.AncillaryCoefficients ancillaryCoefficientsButane(
    n_vapor={-7.17616903,2.53635336,-2.07532869,-2.82241113},
    theta_vapor={1.0,1.5,2.0,4.5},
    n_dliq={1.97874515,0.856799510,-0.341871887,0.304337558},
    theta_dliq={0.345,1.0,1.5,3.0},
    n_dvap={-2.07770057,-3.08362490,-0.485645266,-3.83167519},
    theta_dvap={0.345,5/6,19/6,25/6})
  "Coefficients for the ancillary equations";

  final constant PartialHelmholtzFluid.ThermalConductivityCoefficients thermalConductivityCoefficientsButane(
    reducingTemperature=425.16,
    reducingDensity=3920*0.0581222,
    lambda_0_coeffs=[
     1.62676E-03,    0;
     9.75703E-04,    1;
     2.89887E-02,    2],
    lambda_r_coeffs=[
    -3.04337E-02,    4.18357E-02,    1;
     1.65820E-01,   -1.47163E-01,    2;
    -1.48144E-01,    1.33542E-01,    3;
     5.25500E-02,   -4.85489E-02,    4;
    -6.29367E-03,    6.44307E-03,    5],
    nu=0.63,
    gamma=1.239,
    R0=1.03,
    z=0.063,
    c=1,
    xi_0=0.194E-9,
    Gamma_0=0.0496,
    qd_inverse=0.875350E-9,
    T_ref=637.68) "Coefficients for the thermal conductivity";

  final constant PartialHelmholtzFluid.DynamicViscosityCoefficients dynamicViscosityCoefficientsButane(
    criticalTemperature=425.125,
    criticalMolarVolume=1/3920,
    molarMass=0.0581222,
    epsilon_kappa=280.51,
    sigma=0.57335,
    CET=[
    0.1628213, 0.5],
    a=[
     0.17067154,    0;
    -0.48879666,    1;
     0.039038856,   2],
    b=[
    -19.572881,       0.0;
     219.73999,      -0.25;
    -1015.3226,      -0.5;
     2471.01251,     -0.75;
    -3375.1717,      -1.0;
     2491.6597,      -1.25;
    -787.26086,      -1.5;
     14.085455,      -2.5;
    -0.34664158,     -5.5],
    g=[
     2.30873963359,      0.0,    0,  0,  0;
     2.03404037254,      0.5,    0,  0,  0],
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

    final constant PartialHelmholtzFluid.SurfaceTensionCoefficients surfaceTensionCoefficientsButane(
    coeffs=[
      0.05418,    1.26]) "Coefficients for the surface tension";


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
