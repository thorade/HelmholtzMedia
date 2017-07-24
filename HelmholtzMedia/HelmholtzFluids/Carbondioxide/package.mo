within HelmholtzMedia.HelmholtzFluids;
package Carbondioxide "Carbondioxide"
extends Interfaces.PartialHelmholtzMedium(
  mediumName="Carbondioxide" "short name",
  fluidConstants={fluidConstantsCarbondioxide},
  fluidLimits=fluidLimitsCarbondioxide,
  helmholtzCoefficients=helmholtzCoefficientsCarbondioxide,
  final thermalConductivityCoefficients=thermalConductivityCoefficientsCarbondioxide,
  final dynamicViscosityCoefficients=dynamicViscosityCoefficientsCarbondioxide,
  final surfaceTensionCoefficients=surfaceTensionCoefficientsCarbondioxide,
  final ancillaryCoefficients=ancillaryCoefficientsCarbondioxide,
  Density(min=fluidLimitsCarbondioxide.DMIN, max=fluidLimitsCarbondioxide.DMAX, start=fluidConstantsCarbondioxide.molarMass/fluidConstantsCarbondioxide.criticalMolarVolume),
  Temperature(min=fluidLimitsCarbondioxide.TMIN, max=fluidLimitsCarbondioxide.TMAX, start=298.15),
  AbsolutePressure(min=0, max=801e6, start=101325),
  SpecificEnthalpy(min=fluidLimitsCarbondioxide.HMIN, max=fluidLimitsCarbondioxide.HMAX, start=(fluidLimitsCarbondioxide.HMIN+fluidLimitsCarbondioxide.HMAX)/2),
  SpecificEntropy(min=fluidLimitsCarbondioxide.SMIN, max=fluidLimitsCarbondioxide.SMAX, start=(fluidLimitsCarbondioxide.SMIN+fluidLimitsCarbondioxide.SMAX)/2));

  final constant FluidConstants
  fluidConstantsCarbondioxide(
    casRegistryNumber="124-38-9",
    iupacName="Carbondioxide" "full name",
    structureFormula="CO2",
    chemicalFormula="CO2",
    molarMass=0.0440098,
    gasConstant=8.31451,
    hasCriticalData=true,
       criticalTemperature=304.1282,
       criticalPressure=7377300,
       criticalMolarVolume=0.0440098/467.6,
       HCRIT0=332245.650784287,
       SCRIT0=1433.62533748809,
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
       acentricFactor=0.22394,
    triplePointTemperature=216.592,
    triplePointPressure=0.51795e6,
    normalBoilingPoint=194.686,
    meltingPoint=216.592) "Fluid Constants";

  final constant FluidLimits
  fluidLimitsCarbondioxide(
    TMIN=fluidConstantsCarbondioxide.triplePointTemperature,
    TMAX=2000,
    DMIN=Modelica.Constants.small,
    DMAX=1638.9,
    PMIN=Modelica.Constants.small,
    PMAX=800e6,
    HMIN=-100e3,
    HMAX=+1300e3,
    SMIN=0.52129,
    SMAX=Modelica.Constants.inf) "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsCarbondioxide(
    idealLog=[
      +2.50000,         1.],
    idealPower=[
      -6.124871062444266,         0;
       5.115596317997826,         1],
    idealEinstein=[
      +1.994270420,       -3.15163;
      +0.621052475,       -6.11190;
      +0.411952928,       -6.77708;
      +1.040289220,      -11.32384;
      +0.0832767753,     -27.08792],
    residualPoly=[
       0.388568232032E+00,  0.000,   1.00;
       0.293854759427E+01,  0.750,   1.00;
      -0.558671885349E+01,  1.000,   1.00;
      -0.767531995925E+00,  2.000,   1.00;
       0.317290055804E+00,  0.750,   2.00;
       0.548033158978E+00,  2.000,   2.00;
       0.122794112203E+00,  0.750,   3.00],
    residualBwr=[
       0.216589615432E+01,  1.500,   1.00,    1;
       0.158417351097E+01,  1.500,   2.00,    1;
      -0.231327054055E+00,  2.500,   4.00,    1;
       0.581169164314E-01,  0.000,   5.00,    1;
      -0.553691372054E+00,  1.500,   5.00,    1;
       0.489466159094E+00,  2.000,   5.00,    1;
      -0.242757398435E-01,  0.000,   6.00,    1;
       0.624947905017E-01,  1.000,   6.00,    1;
      -0.121758602252E+00,  2.000,   6.00,    1;
      -0.370556852701E+00,  3.000,   1.00,    2;
      -0.167758797004E-01,  6.000,   1.00,    2;
      -0.119607366380E+00,  3.000,   4.00,    2;
      -0.456193625088E-01,  6.000,   4.00,    2;
       0.356127892703E-01,  8.000,   4.00,    2;
      -0.744277271321E-02,  6.000,   7.00,    2;
      -0.173957049024E-02,  0.000,   8.00,    2;
      -0.218101212895E-01,  7.000,   2.00,    3;
       0.243321665592E-01, 12.000,   3.00,    3;
      -0.374401334235E-01, 16.000,   3.00,    3;
       0.143387157569E+00, 22.000,   5.00,    4;
      -0.134919690833E+00, 24.000,   5.00,    4;
      -0.231512250535E-01, 16.000,   6.00,    4;
       0.123631254929E-01, 24.000,   7.00,    4;
       0.210583219729E-02,  8.000,   8.00,    4;
      -0.339585190264E-03,  2.000,  10.00,    4;
       0.559936517716E-02, 28.000,   4.00,    5;
      -0.303351180556E-03, 14.000,   8.00,    6],
     residualGauss=[
      -0.213654886883E+03,  1.000,   2.00,    2, 2,  -25.,  -325.,  1.16, 1.;
       0.266415691493E+05,  0.000,   2.00,    2, 2,  -25.,  -300.,  1.19, 1.;
      -0.240272122046E+05,  1.000,   2.00,    2, 2,  -25.,  -300.,  1.19, 1.;
      -0.283416034240E+03,  3.000,   3.00,    2, 2,  -15.,  -275.,  1.25, 1.;
       0.212472844002E+03,  3.000,   3.00,    2, 2,  -20.,  -275.,  1.22, 1.],
     residualNonAnalytical=[
      -0.666422765408E+00,  0.000,   1.00,    2, 2,  0.875,  0.300, 0.70, 10.0, 275., 0.3, 3.5;
       0.726086323499E+00,  0.000,   1.00,    2, 2,  0.925,  0.300, 0.70, 10.0, 275., 0.3, 3.5;
       0.550686686128E-01,  0.000,   1.00,    2, 2,  0.875,  0.300, 0.70, 12.5, 275., 1.0, 3.])
  "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsCarbondioxide(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=251.196,
    reducingThermalConductivity_0=0.001,
    lambda_0_num_coeffs=[
       7.5378307E+0,   0.5;
       4.8109652E-2, -99.0],
    lambda_0_den_coeffs=[
       0.4226159E+0,   0.0;
       0.6280115E+0,  -1.0;
      -0.5387661E+0,  -2.0;
       0.6735941E+0,  -3.0;
      -0.4362677E+0,  -6.0;
       0.2255388E+0,  -7.0],
    reducingTemperature_background=1.0,
    reducingMolarVolume_background=1/2272.221E-2,
    reducingThermalConductivity_background=1.0E-3,
    lambda_b_coeffs=[
     2.447164E-2,    0.0,   1.0,   0;
     8.705605E-05,   0.0,   2.0,   0;
    -6.547950E-08,   0.0,   3.0,   0;
     6.594919E-11,   0.0,   4.0,   0],
    gamma=1.2415,
    R0=1.01,
    z=0.065,
    xi_0=1.5E-10,
    Gamma_0=0.052,
    qd_inverse=0.40E-9,
    T_ref=450.0) "Coefficients for the thermal conductivity";

final constant Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsCarbondioxide(
  dynamicViscosityModel=DynamicViscosityModel.VS1,
  collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=1.0,
    epsilon_kappa=251.196,
    CET=[
     1.00697, 0.5],
    a=[
     0.235156,     0;
    -0.491266,     1;
     5.211155E-2,  2;
     5.347906E-2,  3;
    -1.537102E-2,  4],
    b=fill(0.0, 0, 2),
    reducingTemperature_residual=251.196,
    reducingMolarVolume_residual=1/22.7222,
    reducingViscosity_residual=1.0,
    g=fill(0.0, 0, 2),
    e=[
     0.4071119E-2,    0.00,  1.00,  0.00,  0;
     0.7198037E-4,    0.00,  2.00,  0.00,  0;
     0.2411697E-16,  -3.00,  6.00,  0.00,  0;
     0.2971072E-22,   0.00,  8.00,  0.00,  0;
    -0.1627888E-22,  -1.00,  8.00,  0.00,  0],
    nu_po=fill(0.0, 0, 5),
    de_po=fill(0.0, 0, 5))
  "Coefficients for the dynamic viscosity";

  final constant Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsCarbondioxide(
    coeffs=[
      0.05418,    1.26]) "Coefficients for the surface tension";

final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsCarbondioxide(
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -7.0602087,   1.0;
       1.9391218,   1.5;
      -1.6463597,   2.0;
      -3.2995634,   4.5],
    densityLiquidModel=DensityLiquidModel.DL4,
    densityLiquid=[
       1.92451080,   1.02;
      -0.62385555,   1.50;
      -0.32731127,   5.00;
       0.39245142,   5.50],
    densityVaporModel=DensityVaporModel.DV4,
    densityVapor=[
      -1.7074879,   1.02;
      -0.8227467,   1.50;
      -4.6008549,   3.00;
      -10.111178,   7.00;
      -29.742252,  14.00],
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=216.592,
    p_reducing=0.51795e6,
    pressureMelting1=[
       1955.5390,   1;
       2055.4593,   2],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2))
  "Coefficients for the ancillary equations (PS5, DL4, DV4, ML1)";

  annotation (Documentation(info="<html>
These are the coefficients for Carbondioxide.

<dl>
<dt> Span, R. and Wagner, W.</dt>
<dd> <b>A new Equation of State for Carbondioxide covering the Fluid Region from the Triple-Point Temperature to 1100K at Pressures up to 800MPa</b><br />
     Journal of Physical and Chemical Reference Data 25.6, 1509-1596 (1996)<br />
     DOI: <a href=\"https://doi.org/10.1063/1.555991\">10.1063/1.555991</a>
</dd>
<dt> Fenghour, A. and Wakeham, W.A. and Vesovic, V.</dt>
<dd> <b>The Viscosity of Carbon Dioxide</b><br />
     Journal of Physical and Chemical Reference Data 27.1, 31-44 (1998)<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.556013\">10.1063/1.556013</a>
</dd>
<dt> Perkins, Richard A. et. al.</dt>
<dd> <b>Measurement and Correlation of the Thermal Conductivity of Carbondioxide from 135 K to 600 K at Pressures to 70 MPa</b><br />
     Journal of Chemical &amp; Engineering Data 47.5, S. 1263-1271. (2002)<br />
     DOI: <a href=\"http://dx.doi.org/10.1021/je0101202\">10.1021/je0101202</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Carbondioxide;
