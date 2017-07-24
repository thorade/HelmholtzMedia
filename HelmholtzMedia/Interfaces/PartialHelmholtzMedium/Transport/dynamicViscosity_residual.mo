within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function dynamicViscosity_residual
  "Returns dynamic Viscosity residual contribution"
  input ThermodynamicState state;
  output DynamicViscosity eta_r;

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
  Real d_gcm=state.d/1000 "density in g/cm3";

  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_residual=MM/dynamicViscosityCoefficients.reducingMolarVolume_residual;
  Real delta=0 "reduced density";
  Real delta_exp=0 "reduced density in exponential term";
  Real delta_0=0 "close packed density";

  Real[size(dynamicViscosityCoefficients.c, 1),1] c=dynamicViscosityCoefficients.c;

  Real[size(dynamicViscosityCoefficients.g, 1),2] g=dynamicViscosityCoefficients.g;
  Real[size(dynamicViscosityCoefficients.e, 1),5] e=dynamicViscosityCoefficients.e;
  Real[size(dynamicViscosityCoefficients.nu_po, 1),5] nu_po=dynamicViscosityCoefficients.nu_po;
  Real[size(dynamicViscosityCoefficients.de_po, 1),5] de_po=dynamicViscosityCoefficients.de_po;
  // Real[size(dynamicViscosityCoefficients.nu_ex,1),5] nu_ex=dynamicViscosityCoefficients.nu_ex;
  // Real[size(dynamicViscosityCoefficients.de_ex,1),5] de_ex=dynamicViscosityCoefficients.de_ex;

  Real visci=0 "RefProp      visci temporary variable";
  Real xnum=0 "RefProp   numerator temporary variable";
  Real xden=0 "RefProp denominator temporary variable";
  Real G=0 "RefProp temporary variable";
  Real H=0 "RefProp temporary variable";
  Real F=0 "RefProp temporary variable";

  Real x=0 "RefProp temporary variable";
  Real B=0 "RefProp temporary variable";
  Real C=0 "RefProp temporary variable";
  Real D=0 "RefProp temporary variable";
  Real eta_e= 0 "RefProp temporary variable";
  Real eta_0= 0 "RefProp temporary variable";
  Real eta_0a=0 "RefProp temporary variable";
  Real eta_0b=0 "RefProp temporary variable";

  Real eta_red_residual=dynamicViscosityCoefficients.reducingViscosity_residual;

algorithm
  // residual contribution

  if (dynamicViscosityModel == DynamicViscosityModel.VS0) then
    // hardcoded models, use mediumName to distinguish further
    if mediumName == "helium" then
      x := if (state.T<=300.0E0) then log(state.T) else 5.703782474656201;
      B := -47.5295259E0/x + 87.6799309E0 - 42.0741589E0*x + 8.33128289E0*x^2 - 0.589252385E0*x^3;
      C :=  547.309267E0/x - 904.870586E0 + 431.404928E0*x - 81.4504854E0*x^2 + 5.37008433E0*x^3;
      D := -1684.39324E0/x + 3331.08630E0 - 1632.19172E0*x + 308.804413E0*x^2 - 20.2936367E0*x^3;
      eta_e  := exp(min(100, B*d_gcm + C*d_gcm^2 + D*d_gcm^3));
      eta_0a := exp(-0.135311743E0/x + 1.00347841E0 + 1.20654649E0*x - 0.149564551E0*x^2 + 0.0125208416E0*x^3);
      if (state.T < 100.0E0) then
        eta_r := eta_0a*eta_e;
      else
        eta_0b := 196.0E0*state.T^0.71938E0*exp(12.451E0/state.T - 295.67E0/state.T^2 - 4.1249E0);
        if (state.T <=110.0E0) then
          //linear interpolation between T=100...110 K
          eta_0 := eta_0a + (eta_0b-eta_0a)*(state.T-100.0E0)/10.0E0;
        else
          // above 110K use eta_0b
          eta_0 := eta_0b;
        end if;
        eta_r :=eta_0a*eta_e + eta_0 - eta_0a;
      end if;
      eta_r := eta_r/10.0E0;
    elseif mediumName == "water" then
      assert(false, "water viscosity not yet implemented");
    end if;

  elseif ((dynamicViscosityModel == DynamicViscosityModel.VS1)
      or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)) then
    // use the reduced close-packed density delta_0,
    // a simple polynominal, a rational polynominal and an exponential term
    tau := state.T/T_red_residual;
    delta := state.d/d_red_residual;
    if (abs(d_red_residual - 1) > 0.001) then
      delta_exp := state.d/d_crit;
    else
      delta_exp := delta;
    end if;

    if size(g, 1) > 0 then
      // generalized RefProp algorithm, be careful with coeffs: they may differ from article
      if (dynamicViscosityModel == DynamicViscosityModel.VS1) then
        delta_0 := sum(g[i, 1]*tau^g[i, 2] for i in 1:size(g, 1));
      elseif (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative) then
        // alternative inverse form
        delta_0 := g[1, 1]/(1 + sum(g[i, 1]*tau^g[i, 2] for i in 2:size(g, 1)));
      end if;
    else
      delta_0 := 1.0;
    end if;
    for i in 1:size(e, 1) loop
      visci := e[i, 1]*tau^e[i, 2]*delta^e[i, 3]*delta_0^e[i, 4]; // simple polynominal terms
      if (e[i, 5] > 0) then
        visci := visci*exp(-delta_exp^e[i, 5]);
      end if;
      eta_r := eta_r + visci;
    end for;

    for i in 1:size(nu_po, 1) loop
      // numerator of rational poly terms, RefProp algorithm
      xnum := xnum + (nu_po[i, 1]*tau^nu_po[i, 2]*delta^nu_po[i, 3]*delta_0^nu_po[i, 4]);
      if (nu_po[i, 5] > 0) then
        xnum := xnum*exp(-delta_exp^nu_po[i, 5]);
      end if;
    end for;
    for i in 1:size(de_po, 1) loop
      // denominator of rational poly terms, RefProp algorithm
      xden := xden + (de_po[i, 1]*tau^de_po[i, 2]*delta^de_po[i, 3]*delta_0^de_po[i, 4]);
      if (de_po[i, 5] > 0) then
        xden := xden*exp(-delta_exp^de_po[i, 5]);
      end if;
    end for;
    if xden <> 0 then
      eta_r := eta_r + xnum/xden;
    end if;
    // exponential terms not yet implemented!!

  elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
    G := c[1,1] + c[2,1]/state.T;
  //H := sqrt(dm)*(dm-dm_crit)/dm_crit;
    H := sqrt(dm)*(dm- c[8,1])/c[8,1];
    F := G + (c[3,1] + c[4,1]*state.T^(-3/2))*dm^0.1
           + (c[5,1] + c[6,1]/state.T + c[7,1]/state.T^2)*H;
    eta_r :=exp(F) - exp(G);

  elseif (dynamicViscosityModel == DynamicViscosityModel.VS4) then
    assert(false, "VS4 not yet implemented!!");
  else
    assert(false, "unknown dynamicViscosityModel");
  end if;
  eta_r := eta_r*eta_red_residual;

  /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        T = " + String(state.T));
  Modelica.Utilities.Streams.print("      tau = " + String(tau));
  Modelica.Utilities.Streams.print("        d = " + String(state.d));
  Modelica.Utilities.Streams.print("       dm = " + String(dm));
  Modelica.Utilities.Streams.print("    d_gcm = " + String(d_gcm));
  Modelica.Utilities.Streams.print("    delta = " + String(delta));
  Modelica.Utilities.Streams.print("delta_exp = " + String(delta_exp));
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("  delta_0 = " + String(delta_0));
  Modelica.Utilities.Streams.print("     xnum = " + String(xnum) + " and xden = " + String(xden));
  Modelica.Utilities.Streams.print("    eta_r = " + String(eta_r));
  Modelica.Utilities.Streams.print("===========================================");
  */

  annotation (Documentation(info="<html>
<p>
This model is identical to the RefProp VS1 or VS2 model.

The viscosity is split into three contributions:
zero density (dilute gas) viscosity eta_0,
initial density contribution eta_1
and residual contribution eta_r.

This allows to develop functions for each contribution seperately.
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
end dynamicViscosity_residual;
