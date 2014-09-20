within HelmholtzMedia.Examples.Validation;
partial model CheckDerivativesMedium
   replaceable package Medium =
      Modelica.Media.Interfaces.PartialTwoPhaseMedium;

   parameter Medium.AbsolutePressure p0, p1;
   parameter Medium.Temperature T0, T1;
   final parameter Medium.SpecificEnthalpy h0 = Medium.specificEnthalpy_pT(p0,T0);
   final parameter Medium.SpecificEnthalpy h1 = Medium.specificEnthalpy_pT(p0,T1);
   final parameter Modelica.SIunits.Time dt=1;
   final parameter Real dpdt = (p1-p0)/dt;
   final parameter Real dhdt = (h1-h0)/dt;

   Medium.AbsolutePressure p = p0 + dpdt*time;
   Medium.SpecificEnthalpy h = h0 + dhdt*time;

   Medium.ThermodynamicState state_p, state_h;

   Medium.Density d_p, d_h, d_p_int, d_h_int;

   Medium.Density d_p_diff = d_p - d_p_int;
   Medium.Density d_h_diff = d_h - d_h_int;
   Modelica.SIunits.PerUnit d_p_err = d_p_diff/d_p;
   Modelica.SIunits.PerUnit d_h_err = d_h_diff/d_h;

   Real ddph, ddhp;
equation

   state_p = Medium.setState_ph(p,h0);
   state_h = Medium.setState_ph(p0,h);
   d_p = Medium.density(state_p);
   d_h = Medium.density(state_h);

   ddph = Medium.density_derp_h(state_p);
   ddhp = Medium.density_derh_p(state_h);

   der(d_p_int) = ddph*dpdt;
   der(d_h_int) = ddhp*dhdt;
initial equation
   d_p_int = d_p;
   d_h_int = d_h;
end CheckDerivativesMedium;
