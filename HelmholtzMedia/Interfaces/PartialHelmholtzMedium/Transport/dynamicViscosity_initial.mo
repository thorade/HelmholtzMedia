within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function dynamicViscosity_initial
  "Returns dynamic Viscosity initial contribution"
  input ThermodynamicState state;
  output DynamicViscosity eta_1;

protected
  DynamicViscosityModel dynamicViscosityModel=dynamicViscosityCoefficients.dynamicViscosityModel;
  CollisionIntegralModel collisionIntegralModel=dynamicViscosityCoefficients.collisionIntegralModel;

  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Temperature T_red_0=dynamicViscosityCoefficients.reducingTemperature_0;
  Temperature T_red_residual=dynamicViscosityCoefficients.reducingTemperature_residual;
  Real T_star "reduced temperature";
  Real tau "reduced temperature";

  MolarMass MM = fluidConstants[1].molarMass;
  Real dm=state.d/(1000*MM) "molar density in mol/l";     // 1 m3=1000 l
  //Real dm_crit=d_crit/(1000*MM) "molar density in mol/l"; // 1 m3=1000 l

  Real[size(dynamicViscosityCoefficients.b, 1),2] b=dynamicViscosityCoefficients.b;
  Real B_star=0 "reduced second viscosity virial coefficient";
  Real B=0 "second viscosity virial coefficient, l/mol";
  Real sigma=dynamicViscosityCoefficients.sigma;

  Real eta_red_1=dynamicViscosityCoefficients.reducingViscosity_1;

algorithm
  // inital density contribution
  // initialize output
  eta_1 := 0;

  if (dynamicViscosityModel == DynamicViscosityModel.VS0) then
    eta_1 := 0;
  elseif ((dynamicViscosityModel == DynamicViscosityModel.VS1)
      or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)
      or  (dynamicViscosityModel == DynamicViscosityModel.VS4)) then
    // use the second viscosity virial coefficient B according to Rainwater and Friend theory
    T_star := (state.T/dynamicViscosityCoefficients.epsilon_kappa);
    B_star := sum(b[i, 1]*T_star^b[i, 2] for i in 1:size(b, 1));
    B := B_star*0.6022137*sigma^3;
    eta_1 := dynamicViscosity_dilute(state)*B*dm;
  elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
    eta_1 := dm * (b[1,1] + b[2,1]*(b[3,1]-log(state.T/b[4,1]))^2);
  else
    assert(false, "unknown dynamicViscosityModel");
  end if;
  eta_1 := eta_1*eta_red_1;

  /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        T = " + String(state.T));
  Modelica.Utilities.Streams.print("   T_star = " + String(T_star));
  Modelica.Utilities.Streams.print("      tau = " + String(tau));
  Modelica.Utilities.Streams.print("        d = " + String(state.d));
  Modelica.Utilities.Streams.print("       dm = " + String(dm));
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("   B_star = " + String(B_star));
  Modelica.Utilities.Streams.print("        B = " + String(B));
  Modelica.Utilities.Streams.print("    eta_1 = " + String(eta_1));
  Modelica.Utilities.Streams.print("===========================================");
  */

  annotation (Documentation(info="<html>
<p>
This model is identical to the RefProp VS1 or VS2 model.

The viscosity is split into three contributions:
zero density (dilute gas) viscosity eta_0,
initial density contribution eta_1
and residual contribution eta_r.

This allows to develop functions for each contribution separately.
The so called background viscosity is the sum of initial and residual viscosity.

At the critical point and a small region around the critical point, the viscosity is enhanced.
As this critical enhancement is small, it is neglected here.

Special thanks go to Eric W. Lemmon for answering all my emails
and programming a special version of RefProp that outputs also intermediate values.

</p>

<dl>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
<dt>Vogel, E.; K&uuml;chenmeister, C. and Birch, E.</dt>
<dd> <b>Reference correlation of the viscosity of propane</b>.<br />
     Journal of Thermophysics (1998) 10, 417-426.<br />
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
end dynamicViscosity_initial;
