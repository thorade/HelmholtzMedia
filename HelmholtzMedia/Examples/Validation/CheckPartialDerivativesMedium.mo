within HelmholtzMedia.Examples.Validation;
model CheckPartialDerivativesMedium
   replaceable package Medium =
      HelmholtzMedia.HelmholtzFluids.Butane;

   // example for calculating partial derivative wrt h at constant p

   parameter Medium.AbsolutePressure p = 1e5;
   constant Medium.SpecificEnthalpy HMIN= 250e3;
   constant Medium.SpecificEnthalpy HMAX = 350e3;
   Modelica.Blocks.Sources.Sine h_sine(
    freqHz=100,
    startTime=0,
    amplitude=(HMAX - HMIN)/2,
    offset=(HMAX - HMIN)/2 + HMIN)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
   Medium.SpecificEnthalpy h = h_sine.y;

   // state and properties from parameter p and variable h
   Medium.ThermodynamicState state = Medium.setState_ph(p,h);
   Medium.Density d = Medium.density_ph(p, h);

    // time derivatives
    Real dt = der(d);
    Real ht = der(h);

    // partial derivative, two ways
    Real ddhp1 = dt/ht;
    Real ddhp2 = Medium.density_derh_p(state);

    // error, for plotting
    Real absdev = ddhp1-ddhp2;

equation

end CheckPartialDerivativesMedium;
