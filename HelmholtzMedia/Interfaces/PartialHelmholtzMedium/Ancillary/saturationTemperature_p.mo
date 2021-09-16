within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_p
  "ancillary iterative function: calculate saturation temperature for a given pressure by iterating the ancillary function"

  input Modelica.Units.SI.AbsolutePressure p;
  output Modelica.Units.SI.Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real tau "inverse reduced temperature";
  Real T_theta;
  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  Integer nPressureSaturation=size(ancillaryCoefficients.pressureSaturation, 1);
  Real[nPressureSaturation] n=ancillaryCoefficients.pressureSaturation[:, 1];
  Real[nPressureSaturation] theta=ancillaryCoefficients.pressureSaturation[:, 2];

  Real RES_p;
  Real dpdT;
  constant Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
  constant Real tolerance=1e-6 "relative tolerance for RES_p";
  Integer iter=0;
  constant Integer iter_max = 200;

algorithm
if (p<p_crit) and (p>p_trip) then
  // calculate start value from the log(p) vs. 1/T diagram
  // see Span (2000) page 52 / equation 3.98
  T := 1/(1/T_crit - (1/T_trip-1/T_crit)/log(p_crit/p_trip)*log(p/p_crit));
  T := min(T,T_crit- Modelica.Constants.eps);

  // calculate RES_p
  tau := T_crit/T;
  T_theta := max((1 - T/T_crit), Modelica.Constants.small);
  RES_p   := p_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation)) - p;

  while ((abs(RES_p/p)>tolerance) and (iter<iter_max)) loop
    iter:=iter + 1;

    // calculate gradient of RES_p (= gradient of Wagner equation)
    dpdT := -p_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation))
            *(1/T*sum(n[i]*theta[i]*T_theta^(theta[i]-1) for i in 1:nPressureSaturation)
            +tau/T*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation));

    /* // print for debugging
    Modelica.Utilities.Streams.print(" ", "printlog.txt");
    Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
    Modelica.Utilities.Streams.print("T=" + String(T) + " and dpdT=" + String(dpdT), "printlog.txt"); */

    // calculate better T
    T := T - gamma/dpdT*RES_p;

    // check bounds
    T := max(T,T_trip*0.99);
    T := min(T,T_crit);

    // calculate new RES_p
    tau := T_crit/T;
    T_theta := max((1 - T/T_crit), Modelica.Constants.small);
    RES_p := p_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation)) - p;
  end while;
  // Modelica.Utilities.Streams.print("Ancillary.saturationTemperature_p total iteration steps " + String(iter), "printlog.txt");
  assert(iter<iter_max, "Ancillary.saturationTemperature_p did not converge, input was p=" + String(p));

elseif (p<=p_trip) then
  T := T_trip;
elseif (p>=p_crit) then
  T := T_crit;
else
  assert(false, "Ancillary.saturationTemperature_p: this should not happen, check p");
end if;

  // this is an iterative backward function
  // the corresponding ancillary forward function is saturationPressure(T)

  annotation (inverse(p=saturationPressure_T(T=T)));
end saturationTemperature_p;
