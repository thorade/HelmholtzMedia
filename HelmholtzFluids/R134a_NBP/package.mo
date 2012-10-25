within HelmholtzMedia.HelmholtzFluids;
package R134a_NBP "R134a with NBP reference state"
  extends HelmholtzMedia.HelmholtzFluids.R134a(
    helmholtzCoefficients=helmholtzCoefficientsR134a_NBP);

  final constant Interfaces.PartialHelmholtzMedium.EoS.HelmholtzCoefficients
  helmholtzCoefficientsR134a_NBP(
    idealLog=[
      -1.629789E+0,     1.00E0],
    idealPower=[
       9.63437581475882623335,     0.00E0;
       3.61631654246604716326,     1.00E0;
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

end R134a_NBP;
