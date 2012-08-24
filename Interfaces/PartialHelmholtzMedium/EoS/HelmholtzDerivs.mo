within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
record HelmholtzDerivs "dimensionless Helmholtz energy and its derivatives"
  extends Modelica.Icons.Record;

  Temperature T;
  Density d;
  SpecificHeatCapacity R "specific gas constant";
  Real delta(unit="1") "reduced density";
  Real tau(unit="1") "inverse reduced temperature";

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
