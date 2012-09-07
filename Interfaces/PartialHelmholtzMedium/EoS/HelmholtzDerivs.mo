within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
record HelmholtzDerivs "dimensionless Helmholtz energy and its derivatives"
  extends Modelica.Icons.Record;

  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;

  Density d;
  Real delta(unit="1")=d/d_crit "reduced density";
  Temperature T;
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";

  Real i(unit="1") "f_i";
  Real it(unit="1") "(df_i/dtau)@delta=const";
  Real itt(unit="1") "(ddf_i/dtautau)@delta=const";
//Real itd(unit="1") "(ddf_i/dtauddelta)";
//Real id(unit="1") "(df_i/ddelta)@tau=const";
//Real idd(unit="1") "(ddf_i/ddeltadelta)@tau=const";

  Real r(unit="1") "f_r";
  Real rt(unit="1") "(df_r/dtau)@delta=const";
  Real rtt(unit="1") "(ddf_r/dtau)@delta=const";
  Real rtd(unit="1") "(ddf_r/dtaudelta)";
//Real rttd(unit="1") "(df_r/dtautaudelta)";
  Real rd(unit="1") "(df_r/ddelta)@tau=const";
  Real rdd(unit="1") "(ddf_r/ddeltadelta)@tau=const";
end HelmholtzDerivs;
