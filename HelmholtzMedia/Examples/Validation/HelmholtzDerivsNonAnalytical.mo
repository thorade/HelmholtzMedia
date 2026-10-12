within HelmholtzMedia.Examples.Validation;
model HelmholtzDerivsNonAnalytical
  "derivatives of the residual Helmholtz energy of carbon dioxide (with non-analytical terms) against high-precision reference values"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide;
  parameter Real tol=1e-10 "relative tolerance, scaled with max(abs(reference), 1)";
  parameter String fileName="nonanalytical.txt";
  Integer nFailed;

protected
  // reduced density delta and inverse reduced temperature tau, several of them close to the critical point (1, 1)
  constant Real points[12,2]={{0.29999999999999999, 1.2}, {0.90000000000000002, 1.01}, {1.1000000000000001, 0.98999999999999999}, {1.02, 1.0029999999999999}, {0.96999999999999997, 0.998}, {1.5, 0.80000000000000004}, {2.3999999999999999, 1.3}, {1, 1.05}, {1.25, 1}, {0.999, 1.0001}, {1.0009999999999999, 0.99950000000000006}, {0.995, 1}};
  constant Real ref[12,10]={
    {-0.55271366167114921, -1.9554351631710649, -5.1405739462354338, -1.3473126561125912, -4.6846565000552518, -9.3419517715749372, 65.01798741926082, -64.720083030823872, -132.5996929105446, -150.65160273248858},
    {-0.92997844779573346, -0.79266737688113353, 0.50933820903533888, -2.3371311185462895, -19.608691933553228, -2.3556395160697066, -2755.7994591434535, -105.42053550950429, 3.1344694974851341, -0.043115850758168831},
    {-1.0271765281155534, -0.66681027888655009, 0.40531705932266254, -2.5096698422676949, -3.4183782367455784, -1.5711918553590207, -87.208363617073161, 4.4846650561365262, 0.76853690545808662, -0.26817474062015412},
    {-1.0038838832779415, -0.72128598279107703, 0.45123958192017644, -2.4607023943481536, -10.164574077862708, -1.5995056500807268, -1135.5241740603831, 4.3735379955420584, 2.9725042496923915, -0.38007351129442873},
    {-0.95547513477611046, -0.73565956028476986, 0.45673752922265115, -2.3346769684463795, -6.2545468436693845, -1.7127073518601987, -707.25952137637125, -3.036842116314125, 1.9996834946584299, -0.32259904748618956},
    {-0.69746865584821038, -0.1878728583749123, 0.47420754794174752, -2.8398402840148962, -1.3880039160809814, -1.7484715282338432, -0.20740864445668872, 0.014575811412012658, -0.2303970692743838, 0.51549038787301493},
    {-3.0119462168263245, -0.38133439878657416, 1.5880593928993039, -5.2105910177378707, -1.1615455428904162, -2.1263996712473778, 1.2906232356742124, -0.56132198373560205, -0.071318494902572027, 2.923658125435423},
    {-1.4762518446853545, -1.4621570741185088, 19.239792740262669, -48.965220004665866, -5320.562017458753, -86.378340273638429, -550247.16684149136, -9928.7043223297587, 2323.9474421331597, 104.42969385017672},
    {-1.1503062024999671, -0.6232832411469712, 0.37568644905539889, -2.769452166673648, -3.5253064096137083, -1.4371367722334192, -194.62378113293488, 15.152789876712555, -0.42527453912054985, -0.14138436834189003},
    {-0.98163940753180279, -0.72603043595870531, 0.45140260923394759, -2.3998990448307227, -10.123616547421028, -1.6576452290575672, 12499.128885873008, -2.695981838488585, 2.2432714846859754, -0.35279219086727615},
    {-0.98165055492751219, -0.72413691553308057, 0.44940517753172277, -2.39733471751331, -8.0286491951232808, -1.6520619373591616, -2266.3191025204578, -0.13107717817827746, 2.082460438478893, -0.35135098642696028},
    {-0.97849240007065741, -0.72767221572534635, 0.45258233359869482, -2.3920555681176978, -64.55949650133789, -1.6666397896038536, -1316528242.9135613, -10872.31963589465, 2.3875199637480229, -0.34947809014112929}};
  constant String names[10]={"r","rd","rdd","rt","rtt","rtd","rttt","rttd","rtdd","rddd"};
  Real res[10];
  Real dev;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    for i in 1:size(points, 1) loop
      (res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10]) :=
        Medium.EoS.f_r_derivs(delta=points[i,1], tau=points[i,2], order=3);
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
    assert(nFailed == 0, "HelmholtzDerivsNonAnalytical: " + String(nFailed) + " derivatives differ, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks f_r_derivs for carbon dioxide (Span and Wagner 1996), whose equation of state contains three non-analytical terms,
against reference values for all derivatives up to third order, at points around the critical point (delta = 1, tau = 1).</p>
<p>The reference values were computed with mpmath (60 significant digits) by numerical differentiation of the closed-form
residual Helmholtz energy with the coefficients of this library. They agree with
<a href=\"https://github.com/usnistgov/teqp\">teqp</a> (automatic differentiation, CoolProp coefficients) to 1e-10 or better,
except very close to the critical point, where the third derivatives become very large (1e-8), and at delta = 1 exactly,
where the mixed derivatives of teqp deviate by about 2e-5, while the reference values are continuous across delta = 1.</p>
</html>"));
end HelmholtzDerivsNonAnalytical;
