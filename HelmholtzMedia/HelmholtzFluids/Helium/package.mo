within HelmholtzMedia.HelmholtzFluids;
package Helium "Helium"
extends Interfaces.PartialHelmholtzMedium(
  mediumName="helium" "short name",
  fluidConstants={fluidConstantsHelium},
  helmholtzCoefficients=helmholtzCoefficientsHelium,
  thermalConductivityCoefficients=thermalConductivityCoefficientsHelium,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsHelium,
  surfaceTensionCoefficients=surfaceTensionCoefficientsHelium,
  ancillaryCoefficients=ancillaryCoefficientsHelium,
  fluidLimits=fluidLimitsHelium,
  Density(min=fluidLimitsHelium.DMIN, max=fluidLimitsHelium.DMAX, start=fluidConstantsHelium.molarMass/fluidConstantsHelium.criticalMolarVolume),
  Temperature(min=fluidLimitsHelium.TMIN, max=fluidLimitsHelium.TMAX, start=298.15),
  AbsolutePressure(min=0, max=fluidLimitsHelium.PMAX, start=101325),
  SpecificEnthalpy(min=fluidLimitsHelium.HMIN, max=fluidLimitsHelium.HMAX, start=0),
  SpecificEntropy(min=fluidLimitsHelium.SMIN, max=fluidLimitsHelium.SMAX, start=0));

  final constant FluidConstants
  fluidConstantsHelium(
    casRegistryNumber="7440-59-7" "CAS number",
    iupacName="helium-4" "full name",
    structureFormula="He",
    chemicalFormula="He",
    molarMass=0.004002602,
    triplePointTemperature=2.1768,
    normalBoilingPoint=4.223806,
    hasCriticalData=true,
       criticalTemperature=5.19530001363577,
       criticalPressure=228322.86660386,
       criticalMolarVolume=0.004002602/69.5849271894682,
       HCRIT0=11542.3930876446,
       SCRIT0=2115.52517228252,
    hasAcentricFactor=true,
       acentricFactor=-0.3835,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.0,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    triplePointPressure=5039.33,
    meltingPoint=2.1768) "Fluid Constants";

  final constant EoSLimits
  fluidLimitsHelium(
    TMIN=fluidConstantsHelium.triplePointTemperature,
    TMAX=2000,
    DMIN=Modelica.Constants.eps,
    DMAX=565.247,
    PMIN=Modelica.Constants.small,
    PMAX=1000e6,
    HMIN=-261e3,
    HMAX=12945e3,
    SMIN=-4060,
    SMAX=105020) "Helmholtz EoS Limits";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsHelium(
    T_reducing=5.1953,
    d_reducing=0.004002602*17.3837e3,
    useLineSearch=true,
    idealLog=[
          1.5,          1],
    idealPower=[
        0.1733486081588815,  0;
        0.4674523961886186,  1],
    idealEinstein=fill(0.0, 0, 2),
    idealCosh=fill(0.0, 0, 2),
    idealSinh=fill(0.0, 0, 2),
    residualPoly=[
          0.015559018,   1.0,     4;
          3.0638932,     0.425,   1;
         -4.2420844,     0.63,    1;
          0.054418088,   0.69,    2;
         -0.18971904,    1.83,    2;
          0.087856262,   0.575,   3],
    residualBwr=[
          2.2833566,     0.925,   1,  1;
         -0.53331595,    1.585,   1,  2;
         -0.53296502,    1.69,    3,  2;
          0.99444915,    1.51,    2,  1;
         -0.30078896,    2.9,     2,  2;
         -1.6432563,     0.8,     1,  1],
     residualGauss=[
          0.8029102,     1.26,    2,  2, 2,   -1.5497,   -0.2471,    3.15,     0.596;
          0.026838669,   3.51,    1,  2, 2,   -9.245,    -0.0983,    2.54505,  0.3423;
          0.04687678,    2.785,   2,  2, 2,   -4.76323,  -0.1556,    1.2513,   0.761;
         -0.14832766,    1.0,     1,  2, 2,   -6.3826,   -2.6782,    1.9416,   0.9747;
          0.03016211,    4.22,    1,  2, 2,   -8.7023,   -2.7077,    0.5984,   0.5868;
         -0.019986041,   0.83,    3,  2, 2,   -0.255,    -0.6621,    2.2282,   0.5627;
          0.14283514,    1.575,   2,  2, 2,   -0.3523,   -0.1775,    1.606,    2.5346;
          0.007418269,   3.447,   2,  2, 2,   -0.1492,   -0.4821,    3.815,    3.6763;
         -0.22989793,    0.73,    3,  2, 2,   -0.05,     -0.3069,    1.61958,  4.5245;
          0.79224829,    1.634,   2,  2, 2,   -0.1668,   -0.1758,    0.6407,   5.039;
         -0.049386338,   6.13,    2,  2, 2,  -42.2358, -1357.6577,   1.076,    0.959])
  "Coefficients of the Helmholtz EoS";

  final constant Transport.ThermalConductivityCoefficients
  thermalConductivityCoefficientsHelium(
    thermalConductivityModel=ThermalConductivityModel.TC0,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK0,
    reducingTemperature_0=10,
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
  dynamicViscosityCoefficientsHelium(
    dynamicViscosityModel=DynamicViscosityModel.VS0,
    collisionIntegralModel=CollisionIntegralModel.CI0,
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
  surfaceTensionCoefficientsHelium(
    coeffs=[
       0.0004656,   1.04;
       0.001889,    2.468;
      -0.002006,    2.661]) "Coefficients for the surface tension";

  final constant Ancillary.AncillaryCoefficients
  ancillaryCoefficientsHelium(
    T_crit=5.1953,
    p_crit=228320,
    d_crit=0.004002602*17.3837e3,
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=10,
    p_reducing=1000e3,
    pressureMelting1=[
      -1.7455837,      0.000000;
       1.6979793,      1.555414],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2),
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -4.0558,   1.0;
      0.69089,   1.5;
      1.1982,   2.0;
      0.79703,   5.0;
      -9.0432,   14.0],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
      2.697,   0.44;
      -4.0115,   0.72;
      5.5835,   1.0;
      -4.7467,   1.5;
      2.7509,   2.5;
      -2.4711,   5.0],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -1.8975,   0.41;
      -6.9133,   1.35;
      13.711,   2.0;
      -14.644,   2.5;
      -23.179,   10.0])
  "Coefficients for the ancillary equations (PS5, DL1, DV3, ML1)";


  annotation (Documentation(info="<html>
These are the coefficients for Helium.
The reference state is set to NBP (saturated liquid at 101325 Pa: h = 0, s = 0), the REFPROP default.
<p>
The equation of state was fitted with the molar gas constant R = 8.314472 J/(mol K),
whereas HelmholtzMedia uses the exact value Modelica.Constants.R = 8.31446261815324 J/(mol K) for all fluids.
The relative difference of 1.1e-6 is far below the uncertainty of the equation of state;
the check values in Tables 3 and 4 of NIST IR 8474 are reproduced within their stated digits
(see Examples.Validation.EoSHelium).
</p>
<p>
The ancillary equations for saturation pressure and densities are those of REFPROP 10.
</p>

<dl>
<dt> Ortiz Vega, D.O., Hall, K.R., Holste, J.C., Harvey, A.H. and Lemmon, E.W.</dt>
<dd> <b>An Equation of State for the Thermodynamic Properties of Helium</b><br />
     NIST Internal Report 8474, National Institute of Standards and Technology, 2023.<br />
     DOI: <a href=\"https://doi.org/10.6028/NIST.IR.8474\">10.6028/NIST.IR.8474</a>
</dd>
<dt> Arp, V.D., McCarty, R.D. and Friend, D.G.</dt>
<dd> <b>Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa</b> (viscosity)<br />
     NIST Technical Note 1334 (revised), 1998.<br />
     URL: <a href=\"https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf\">nbstechnicalnote1334.pdf</a>
</dd>
<dt> Hands, B.A. and Arp, V.D.</dt>
<dd> <b>A Correlation of Thermal Conductivity Data for Helium</b><br />
     Cryogenics, 21(12):697-703, 1981.<br />
     DOI: <a href=\"https://doi.org/10.1016/0011-2275(81)90211-3\">10.1016/0011-2275(81)90211-3</a>
</dd>
<dt> Mulero, A., Cachadi&ntilde;a, I. and Parra, M.I.</dt>
<dd> <b>Recommended Correlations for the Surface Tension of Common Fluids</b><br />
     J. Phys. Chem. Ref. Data, 41:043105, 2012.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.4768782\">10.1063/1.4768782</a>
</dd>
<dt> McCarty, R.D. and Arp, V.D.</dt>
<dd> <b>A new wide range equation of state for helium</b> (melting line)<br />
     Adv. Cryo. Eng. 35:1465-1475 (1990)
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Helium;
