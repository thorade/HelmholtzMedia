within HelmholtzMedia.Validation.Derivatives;
model HelmholtzDerivsAssociating
  "derivatives of the residual Helmholtz energy of ammonia (with associating terms) against high-precision reference values"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia;
  parameter Real tol=1e-10 "relative tolerance, scaled with max(abs(reference), 1)";
  parameter String fileName="associating.txt";
  Integer nFailed;

protected
  // reduced density delta and inverse reduced temperature tau, including the centres gamma of the associating terms (1.108, 1.313)
  constant Real points[10,2]={{0.3, 1.2}, {0.9, 1.01}, {1.1, 0.99}, {1.0, 1.0}, {1.5, 0.8}, {2.4, 1.3}, {3.0, 2.0}, {0.05, 0.6}, {1.25, 1.108}, {0.45, 1.313}};
  constant Real ref[10,10]={
    {-0.68116178756638213, -2.0746963837879604, 1.4591255195151607, -1.7150140110292382, -4.9320861535167939, -5.6188374670478618, -22.2681147487339, -20.416264434814096, 2.758900946358279, 1.5037799713128197},
    {-1.0432555993081523, -0.82687321438075034, 0.60043664541807251, -2.6771029157159822, -20.145216737387106, -2.1580520057104446, -1531.0326108748438, -57.185784752426848, 5.971702418479965, -0.27091123035186487},
    {-1.1400836168543267, -0.6917261957341817, 0.44625910976268157, -2.7897103506024696, -4.2188601876502103, -1.4542143734599207, -179.38934381778432, 6.3107121582610111, 0.94233424323404635, -0.38469662651126629},
    {-1.0952929419598878, -0.75394904380842742, 0.50789808763785627, -2.7053771301977085, -9.5760327225563167, -1.5277957158478097, -815.4916303021123, 6.9880099771937099, 2.7371525333328903, -0.52369426287140475},
    {-0.78956262124937715, -0.25457399554579679, 0.37782514747819534, -2.9381386510891825, -1.9131339465792412, -1.5167728837777767, -0.81694380016329349, 0.082853102364579115, 0.077537185622398322, 0.2281304011742086},
    {-3.1905272742181987, -0.48554017378368725, 0.90357899798555708, -5.1484594371195573, -1.4010098407302631, -1.5050150411023717, 0.9056782078148351, 0.37787141316613992, -0.43398373226012371, 1.6443839863253387},
    {-7.7751012008597584, -0.61842832226201427, 2.5013971431333366, -6.7870484630271743, -0.73931040317082342, -1.467793413766755, 0.87308009246642337, -0.1499972195685467, 0.52835157489218842, 3.6941965878230265},
    {-0.016603245053784112, -0.32639567176671253, 0.22713160305971786, -0.10338074441173203, -0.15521666658717995, -2.0523385817037103, -0.37641588477363758, -3.0311292844965205, 0.6099064165374809, 0.026631174312176955},
    {-1.6195251169597422, -0.74248272217335032, 0.25891092929971349, -3.3102088050936249, -13.750030068850331, -2.6738706390966707, 112.88191286355038, 115.38336569225686, 16.26913124911016, -0.47316955433460089},
    {-1.31981527520512, -2.5420864477760601, 2.6231122371845474, -3.6843962541911108, -13.068071986787087, -7.7520915298703521, -52.134075081512343, -32.085126602819088, 11.120679416472546, 5.4757414610415411}};
  constant String names[10]={"r","rd","rdd","rt","rtt","rtd","rttt","rttd","rtdd","rddd"};
  Real res[10];
  Real dev;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    for i in 1:size(points, 1) loop
      (res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10]) :=
        Medium.Internal.EoS.f_r_derivs(delta=points[i,1], tau=points[i,2], order=3);
      for k in 1:10 loop
        dev := abs(res[k] - ref[i,k])/max(abs(ref[i,k]), 1);
        if dev > tol then
          nFailed := nFailed + 1;
          Modelica.Utilities.Streams.print("delta=" + String(points[i,1]) + " tau=" + String(points[i,2]) + " " + names[k]
            + " = " + String(res[k], significantDigits=17) + ", reference " + String(ref[i,k], significantDigits=17)
            + ", rel. deviation " + String(dev, significantDigits=3), fileName);
        end if;
      end for;
    end for;
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "HelmholtzDerivsAssociating: " + String(nFailed) + " derivatives differ, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks f_r_derivs for ammonia (Gao et al. 2023), whose equation of state contains two associating terms
n*delta^d*tau^t*exp(eta*(delta-epsilon)^2 + 1/(beta*(tau-gamma)^2 + b)),
against reference values for all derivatives up to third order, including points at the centres gamma of the associating terms.</p>
<p>The reference values were computed with sympy (symbolic differentiation of the closed-form residual Helmholtz energy
with the coefficients of this library, evaluated with 25 significant digits).</p>
</html>"));
end HelmholtzDerivsAssociating;
