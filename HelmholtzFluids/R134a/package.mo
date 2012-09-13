within HelmholtzMedia.HelmholtzFluids;
package R134a "R134a"
  extends Interfaces.PartialHelmholtzMedium(
    fluidConstants={fluidConstantsR134a},
    helmholtzCoefficients=helmholtzCoefficientsR134a,
    thermalConductivityCoefficients=thermalConductivityCoefficientsR134a,
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsR134a,
    surfaceTensionCoefficients=surfaceTensionCoefficientsR134a,
    ancillaryCoefficients=ancillaryCoefficientsR134a,
    fluidLimits=fluidLimitsR134a,
    Density(min=fluidLimitsR134a.DMIN, max=fluidLimitsR134a.DMAX, start=fluidConstantsR134a.molarMass/fluidConstantsR134a.criticalMolarVolume),
    Temperature(min=fluidLimitsR134a.TMIN, max=fluidLimitsR134a.TMAX, start=298.15),
    AbsolutePressure(min=0, max=70e6, start=101325),
    SpecificEnthalpy(min=fluidLimitsR134a.HMIN, max=fluidLimitsR134a.HMAX, start=(fluidLimitsR134a.HMIN+fluidLimitsR134a.HMAX)/2),
    SpecificEntropy(min=fluidLimitsR134a.SMIN, max=fluidLimitsR134a.SMAX, start=(fluidLimitsR134a.SMIN+fluidLimitsR134a.SMAX)/2));

  final constant Interfaces.PartialHelmholtzMedium.FluidConstants
  fluidConstantsR134a(
    chemicalFormula="C2H2F4",
    structureFormula="",
    casRegistryNumber="811-97-2",
    iupacName="",
    molarMass=0.102032,
    hasCriticalData=true,
       criticalTemperature=374.21,
       criticalPressure=4059280,
       criticalMolarVolume=0.102032/511.899952,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=2.058,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.32684,
    triplePointTemperature=169.85,
    triplePointPressure=389.517537,
    normalBoilingPoint=247.076,
    meltingPoint=247.076063) "Fluid Constants";

  final constant Interfaces.PartialHelmholtzMedium.FluidLimits
  fluidLimitsR134a(
    TMIN=fluidConstantsR134a.triplePointTemperature,
    TMAX=455,
    DMIN=Modelica.Constants.small,
    DMAX=1592,
    PMIN=Modelica.Constants.small,
    PMAX=70e6,
    HMIN=0,
    HMAX=600e3,
    SMIN=0,
    SMAX=1e12) "Helmholtz EoS Limits";

  final constant Interfaces.PartialHelmholtzMedium.EoS.HelmholtzCoefficients
  helmholtzCoefficientsR134a(
    idealLog=[
      -1.629789E+0,     1.00E0],
    idealPower=[
      -1.019535E+0,     0.00E0;
       9.047135E+0,     1.00E0;
      -9.723916E+0,    -0.50E0;
      -3.927170E+0,    -0.75E0],
    idealEinstein=fill(0.0, 0, 2),
    residualPoly=[
      0.5586817000E-01,  -0.50,   2.00,   0;
      0.4982230000E+00,   0.00,   1.00,   0;
      0.2458698000E-01,   0.00,   3.00,   0;
      0.8570145000E-03,   0.00,   6.00,   0;
      0.4788584000E-03,   1.50,   6.00,   0;
     -0.1800808000E+01,   1.50,   1.00,   0;
      0.2671641000E+00,   2.00,   1.00,   0;
     -0.4781652000E-01,   2.00,   2.00,   0],
    residualBwr=[
      0.1423987000E-01,   1.00,   5.00,   1;
      0.3324062000E+00,   3.00,   2.00,   1;
     -0.7485907000E-02,   5.00,   2.00,   1;
      0.1017263000E-03,   1.00,   4.00,   2;
     -0.5184567000E+00,   5.00,   1.00,   2;
     -0.8692288000E-01,   5.00,   4.00,   2;
      0.2057144000E+00,   6.00,   1.00,   2;
     -0.5000457000E-02,  10.00,   2.00,   2;
      0.4603262000E-03,  10.00,   4.00,   2;
     -0.3497836000E-02,  10.00,   1.00,   3;
      0.6995038000E-02,  18.00,   5.00,   3;
     -0.1452184000E-01,  22.00,   3.00,   3;
     -0.1285458000E-03,  50.00,  10.00,   4],
   residualGauss=fill(0.0, 0, 9)) "Coefficients of the Helmholtz EoS";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsR134a(
    reducingTemperature_0=1,
    reducingThermalConductivity_0=1,
    lambda_0_num_coeffs=[
    -1.05248E-2,    0;
     8.00982E-5,    1],
    reducingTemperature_residual=1,
    reducingMolarVolume_residual=1/5049.886,
    reducingThermalConductivity_residual=2.055E-3,
    lambda_r_coeffs=[
     1.836526E+0,   0,   1,   0;
     5.126143E+0,   0,   2,   0;
    -1.436883E+0,   0,   3,   0;
     6.261441E-1,   0,   4,   0],
    xi_0=1.94E-10,
    Gamma_0=0.0496,
    qd_inverse=5.285356E-10,
    T_ref=561.411) "Coefficients for the thermal conductivity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsR134a(
    dynamicViscosityModel=Interfaces.PartialHelmholtzMedium.Types.DynamicViscosityModel.VS1_alternative,
    collisionIntegralModel=Interfaces.PartialHelmholtzMedium.Types.CollisionIntegralModel.CI1,
    sigma=0.468932,
    epsilon_kappa=299.363,
    CET=[
     0.215729E0, 0.5],
    a=[
      0.355404E+0,  0;
     -0.464337E+0,  1;
      0.257353E-1,  2],
    b=[
    -0.19572881E+2,   0.00;
     0.21973999E+3,  -0.25;
    -0.10153226E+4,  -0.50;
     0.24710125E+4,  -0.75;
    -0.33751717E+4,  -1.00;
     0.24916597E+4,  -1.25;
    -0.78726086E+3,  -1.50;
     0.14085455E+2,  -2.50;
    -0.34664158E+0,  -5.50],
    reducingTemperature_residual=374.21,
    reducingMolarVolume_residual=1/5017.0613,
    reducingViscosity_residual=1e3,
    g=[
     3.163695635587490,      0.00;
    -0.8901733752064137E-1,  1.00;
     0.1000352946668359,     2.00],
    e=[
    -0.2069007192080741E-1,  0.00,  1.00,  0.00,  0;
     0.3560295489828222E-3, -6.00,  2.00,  0.00,  0;
     0.2111018162451597E-2, -2.00,  2.00,  0.00,  0;
     0.1396014148308975E-1, -0.50,  2.00,  0.00,  0;
    -0.4564350196734897E-2,  2.00,  2.00,  0.00,  0;
    -0.3515932745836890E-2,  0.00,  3.00,  0.00,  0;
    -0.2147633195397038,     0.00,  0.00, -1.00,  0],
    nu_po=[
     0.2147633195397038,     0.00,  0.00,  0.00,  0],
    de_po=[
     1.00,                   0.00,  0.00,  1.00,  0;
    -1.00,                   0.00,  1.00,  0.00,  0])
  "Coefficients for the dynamic viscosity";

  final constant
  Interfaces.PartialHelmholtzMedium.Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsR134a(
    coeffs=[
      0.06016,    1.26]) "Coefficients for the surface tension";

  final constant
  Interfaces.PartialHelmholtzMedium.Ancillary.AncillaryCoefficients
  ancillaryCoefficientsR134a(
    pressureSaturationModel=Interfaces.PartialHelmholtzMedium.Types.PressureSaturationModel.PS5,
    pressureSaturation=[
      -0.77513E+01,   1.0;
       0.29263E+01,   1.5;
      -0.26622E+01,   1.9;
      -0.39711E+01,   4.25],
    densityLiquidModel=Interfaces.PartialHelmholtzMedium.Types.DensityLiquidModel.DL1,
    densityLiquid=[
       0.12449E+02,   0.5;
      -0.41023E+02,   0.7;
       0.73641E+02,   0.9;
      -0.64635E+02,   1.1;
       0.22551E+02,   1.3],
    densityVaporModel=Interfaces.PartialHelmholtzMedium.Types.DensityVaporModel.DV3,
    densityVapor=[
      -0.29174E+01,   0.383;
      -0.72542E+01,   1.21;
      -0.23306E+02,   3.3;
       0.59840E+01,   5.6;
      -0.71821E+02,   7.0]) "Coefficients for the ancillary equations";


  annotation (Documentation(info="<html>
These are the coefficients for R134a. 

<dl>
<dt> Tillner-Roth, R. and Baehr, H.D.,</dt>
<dd> <b>An international standard formulation of the thermodynamic properties of 1,1,1,2-tetrafluoroethane (HFC-134a) for temperatures from 170 K to 455 K at pressures up to 70 MPa</b><br>
     J. Phys. Chem. Ref. Data, 23:657-729, 1994.<br>
     DOI: <a href=\"http://dx.doi.org/10.1063/1.555958\">10.1063/1.555958</a>
</dd>
<dt> Huber, M.L., Laesecke, A., and Perkins, R.A.</dt>
<dd> <b>Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a</b><br>
     Ind. Eng. Chem. Res., 42:3163-3178, 2003.<br>
     DOI: <a href=\"http://dx.doi.org/10.1021/ie0300880\">10.1021/ie0300880</a>
</dd>
<dt> Perkins, R.A., Laesecke, A., Howley, J., Ramires, M.L.V., Gurova, A.N., and Cusco, L.</dt>
<dd> <b>Experimental thermal conductivity values for the IUPAC round-robin sample of 1,1,1,2-tetrafluoroethane (R134a)</b><br>
     NISTIR, 2000.
</dd><dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br>
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br>
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end R134a;
