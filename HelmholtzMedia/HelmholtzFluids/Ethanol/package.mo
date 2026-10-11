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
  AbsolutePressure(min=0, max=fluidLimitsEthanol.PMAX, start=101325),
  SpecificEnthalpy(min=fluidLimitsEthanol.HMIN, max=fluidLimitsEthanol.HMAX, start=(fluidLimitsEthanol.HMIN+fluidLimitsEthanol.HMAX)/2),
  SpecificEntropy(min=fluidLimitsEthanol.SMIN, max=fluidLimitsEthanol.SMAX, start=(fluidLimitsEthanol.SMIN+fluidLimitsEthanol.SMAX)/2));

  final constant FluidConstants
  fluidConstantsEthanol(
    chemicalFormula="C2H6O",
    structureFormula="",
    casRegistryNumber="64-17-5",
    iupacName="",
    molarMass=0.04606844,
    hasCriticalData=true,
       criticalTemperature=514.709284881,
       criticalPressure=6267907.51014,
       criticalMolarVolume=0.04606844/273.217468213,
       HCRIT0=758708.687380538,
       SCRIT0=1680.2329679611,
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
       acentricFactor=0.646,
    triplePointTemperature=159,
    triplePointPressure=7.184e-4,
    normalBoilingPoint=351.57,
    meltingPoint=159) "Fluid Constants";

  final constant EoSLimits
  fluidLimitsEthanol(
    TMIN=fluidConstantsEthanol.triplePointTemperature,
    TMAX=650,
    DMIN=Modelica.Constants.eps,
    DMAX=1031,
    PMIN=Modelica.Constants.small,
    PMAX=280e6,
    HMIN=-442e3,
    HMAX=1487e3,
    SMIN=-2180,
    SMAX=10170) "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsEthanol(
    T_reducing=514.71,
    d_reducing=0.04606844*5930,
    idealLog=[
      3.43069,          1],
    idealPower=[
     -3.644399430821064,          0;
      5.07080732963899,          1],
    idealEinstein=[
      2.14326,  -420.4/514.71;
      5.09206,  -1334/514.71;
      6.60138,  -1958/514.71;
      5.70777,  -4420/514.71],
    residualPoly=[
      0.058200796,    1.00,   4;
      0.94391227,     1.04,   1;
     -0.80941908,     2.72,   1;
      0.55359038,     1.174,  2;
     -1.4269032,      1.329,  2;
      0.13448717,     0.195,  3],
    residualBwr=[
      0.42671978,     2.43,   1,   1;
     -1.1700261,      1.274,  1,   1;
     -0.92405872,     4.16,   1,   2;
      0.34891808,     3.30,   3,   1;
     -0.91327720,     4.177,  3,   2;
      0.022629481,    2.50,   2,   1;
     -0.15513423,     0.81,   2,   2;
      0.21055146,     2.02,   6,   1;
     -0.21997690,     1.606,  6,   1;
     -0.0065857238,   0.86,   8,   1],
    residualGauss=[
      0.75564749,     2.50,   1,   2,   2,   -1.075,   -1.207,   1.194,   0.779;
      0.10694110,     3.72,   1,   2,   2,   -0.463,   -0.0895,  1.986,   0.805;
     -0.069533844,    1.19,   2,   2,   2,   -0.876,   -0.581,   1.583,   1.869;
     -0.24947395,     3.25,   3,   2,   2,   -1.108,   -0.947,   0.756,   0.694;
      0.027177891,    3.00,   3,   2,   2,   -0.741,   -2.356,   0.495,   1.312;
     -0.0009053953,   2.00,   2,   2,   2,   -4.032,   -27.01,   1.002,   2.054;
     -0.12310953,     2.00,   2,   2,   2,   -2.453,   -4.542,   1.077,   0.441;
     -0.08977971,     1.00,   2,   2,   2,   -2.300,   -1.287,   1.493,   0.793;
     -0.39512601,     1.00,   1,   2,   2,   -3.143,   -3.090,   1.542,   0.313]) "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsEthanol(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=514.71,
    reducingThermalConductivity_0=1e-3,
    lambda_0_num_coeffs=[
     -2.09575,     0;
      19.9045,     1;
     -53.964,      2;
      82.1223,     3;
     -1.98864,     4;
     -0.495513,    5],
    lambda_0_den_coeffs=[
      0.172230,    0;
     -0.078273,    1;
      1.0,         2],
    reducingTemperature_background=514.71,
    reducingMolarVolume_background=1/5930,
    reducingThermalConductivity_background=1,
    lambda_b_coeffs=[
      0.0267222,   0,   1,   0;
      0.148279,    0,   2,   0;
     -0.130429,    0,   3,   0;
      0.0346232,   0,   4,   0;
     -0.00244293,  0,   5,   0;
      0.0177166,   1,   1,   0;
     -0.0893088,   1,   2,   0;
      0.0684664,   1,   3,   0;
     -0.0145702,   1,   4,   0;
      0.000809189, 1,   5,   0],
    R0=1.02,
    xi_0=0.164296E-9,
    Gamma_0=0.05885,
    qd_inverse=0.53E-9,
    T_ref=772.06) "Coefficients for the thermal conductivity";

  final constant Transport.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsEthanol(
    dynamicViscosityModel=DynamicViscosityModel.VS1,
    collisionIntegralModel=CollisionIntegralModel.CI1,
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

  final constant Transport.SurfaceTensionCoefficients
  surfaceTensionCoefficientsEthanol(
    coeffs=[
      0.0628594,   1.04726;
     -29.4208,     4.34139;
      29.4913,     4.34996]) "Coefficients for the surface tension";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsEthanol(
    T_crit=514.71,
    p_crit=6268000,
    d_crit=0.04606844*5930,
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
     -9.1043,     1.0;
      4.7263,     1.5;
     -9.7145,     2.0;
      4.1536,     2.55;
     -2.0739,     4.0],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
      11.632,     0.66;
     -218.66,     1.5;
      826.94,     1.9;
     -1351.2,     2.3;
      1051.7,     2.7;
     -318.09,     3.1],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      160.386,    0.49;
     -302.488,    0.55;
      504.347,    0.84;
     -735.211,    1.1;
      868.411,    1.5;
     -528.497,    1.7;
     -55.7818,    6.0;
     -124.911,    13.0])
  "Coefficients for the ancillary equations (PS5, DL1, DV3)";


  annotation (Documentation(info="<html>
These are the coefficients for Ethanol.
The reference state is set to NBP (saturated liquid at 101325 Pa: h = 0, s = 0), the REFPROP default.

<dl>
<dt> Schroeder, J. A., Penoncello, S. G., and Schroeder, J. S.</dt>
<dd> <b>A Fundamental Equation of State for Ethanol</b><br />
     J. Phys. Chem. Ref. Data 43, 043102 (2014)<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.4895394\">10.1063/1.4895394</a>
</dd>
<dt> Kiselev, S. B., Ely, J. F., Abdulagatov, I. M., and Huber, M. L.</dt>
<dd> <b>Generalized SAFT-DFT/DMT Model for the Thermodynamic, Interfacial, and Transport Properties of Associating Fluids: Application for n-Alkanols</b><br />
     Ind. Eng. Chem. Res. 44, 6916-6927 (2005)<br />
     DOI: <a href=\"http://dx.doi.org/10.1021/ie050010e\">10.1021/ie050010e</a>
</dd>
<dt> Assael, M. J., Sykioti, E. A., Huber, M. L., and Perkins, R. A.</dt>
<dd> <b>Reference Correlation of the Thermal Conductivity of Ethanol from the Triple Point to 600 K and up to 245 MPa</b><br />
     J. Phys. Chem. Ref. Data 42, 023102 (2013)<br />
     DOI: <a href=\"http://dx.doi.org/10.1063/1.4797368\">10.1063/1.4797368</a>
</dd>
<dt> Mulero, A., Cachadi&ntilde;a, I., and Sanjuan, E. L.</dt>
<dd> <b>Surface Tension of Alcohols. Data Selection and Recommended Correlations</b><br />
     J. Phys. Chem. Ref. Data 44, 033104 (2015)
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Ethanol;
