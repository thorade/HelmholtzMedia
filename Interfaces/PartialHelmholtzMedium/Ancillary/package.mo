within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
package Ancillary 


  function saturationPressure_T
  "ancillary function: calculate saturation pressure for a given Temperature"
    input Modelica.SIunits.Temperature T;
    output Modelica.SIunits.AbsolutePressure p;

protected
    constant Temperature T_crit=fluidConstants[1].criticalTemperature;
    constant Real tau=T_crit/T "inverse reduced temperature";
    constant Real T_theta=max((1 - T/T_crit), Modelica.Constants.small);
    constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

    Integer nPressureSaturation=size(ancillaryCoefficients.pressureSaturation, 1);
    Real[nPressureSaturation] n=ancillaryCoefficients.pressureSaturation[:, 1];
    Real[nPressureSaturation] theta=ancillaryCoefficients.pressureSaturation[:, 2];

  algorithm
    assert(T <= T_crit, "saturationPressure error: Temperature is higher than critical temperature");
    p := p_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation));

    // this is an ancillary forward function
    // the corresponding iterative backward function is saturationTemperature(p)
    annotation (inverse(T=saturationTemperature_p(p=p)), Documentation(info="<html>
      <p>
      This algorithm returns the saturation pressure as a function of Temperature: psat=psat(T).
      This type of vapor pressure equation was developed by W. Wagner.
      Because it cannot be solved for temperature analytically, 
      the inverse function Tsat=Tsat(p) has to find Tsat iteratively.
      </p>
      
      <dl>
      <dt>Wagner, W.</dt>
      <dd> <b>Eine mathematisch statistische Methode zum Aufstellen thermodynamischer Gleichungen - gezeigt am Beispiel der Dampfdruckkurve reiner fluider Stoffe.</b><br>
           Forschrittberichte der VDI Zeitschriften, Reihe 3, Nr. 39 (1974)
      </dd>
      </dl>
      </html>"));
  end saturationPressure_T;


  function saturationTemperature_p
  "ancillary iterative function: calculate saturation temperature for a given pressure by iterating the ancillary function"

    input Modelica.SIunits.AbsolutePressure p;
    output Modelica.SIunits.Temperature T;

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
    constant Real tolerance=1e-9 "relative tolerance for RES_p";
    Integer iter=0;
    constant Integer iter_max = 200;

  algorithm
    assert(p >= p_trip, "saturationTemperature error: Pressure is lower than triple-point pressure");
    assert(p <= p_crit, "saturationTemperature error: Pressure is higher than critical pressure");

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
      T := max(T,T_trip);
      T := min(T,T_crit);

      // calculate new RES_p
      tau := T_crit/T;
      T_theta := max((1 - T/T_crit), Modelica.Constants.small);
      RES_p := p_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nPressureSaturation)) - p;
    end while;
    // Modelica.Utilities.Streams.print("saturationTemperature_p total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "saturationTemperature_p did not converge, input was p=" + String(p));

    // this is an iterative backward function
    // the corresponding ancillary forward function is saturationPressure(T)
    annotation (inverse(p=saturationPressure_T(T=T)));
  end saturationTemperature_p;

end Ancillary;
