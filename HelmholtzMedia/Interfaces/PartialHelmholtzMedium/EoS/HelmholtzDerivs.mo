within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
record HelmholtzDerivs "dimensionless Helmholtz energy and its derivatives"
  extends Modelica.Icons.R_secord;

  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=fluidConstants[1].gasConstant/MM "specific gas constant";
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;

  Density d(min=0)=d_crit;
  Real delta(unit="1")=d/d_crit "reduced density";
  Temperature T=T_crit;
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";

  Real i(unit="1")=0 "f_i";
//Real id(unit="1")=0 "(df_i/ddelta)@tau=const";
//Real idd(unit="1")=0 "(d2f_i/ddelta2)@tau=const";
  Real it(unit="1")=0 "(df_i/dtau)@delta=const";
  Real itt(unit="1")=0 "(d2f_i/dtau2)@delta=const";
//Real itd(unit="1")=0 "(d2f_i/dtau ddelta)";
  Real ittt(unit="1")=0 "(d3f_i/dtau3)@delta=const";

  Real r(unit="1")=0 "f_r";
  Real rd(unit="1")=0 "(df_r/ddelta)@tau=const";
  Real rdd(unit="1")=0 "(d2f_r/ddelta2)@tau=const";
  Real rt(unit="1")=0 "(df_r/dtau)@delta=const";
  Real rtt(unit="1")=0 "(d2f_r/dtau2)@delta=const";
  Real rtd(unit="1")=0 "(d2f_r/dtau ddelta)";

  Real rttt(unit="1")=0 "(d3f_r/dtau3)@delta=const";
  Real rttd(unit="1")=0 "(d3f_r/dtau2 ddelta)";
  Real rtdd(unit="1")=0 "(d3f_r/dtau ddelta2)";
  Real rddd(unit="1")=0 "(d3f_r/ddelta3)@tau=const";
end HelmholtzDerivs;
