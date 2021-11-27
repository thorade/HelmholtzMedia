within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationPressure_T
  "ancillary function: calculate saturation pressure for a given Temperature"
  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.AbsolutePressure p;

protected
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real tau=T_crit/T "inverse reduced temperature";
  Real T_theta=max((1 - T/T_crit), Modelica.Constants.small);
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  Integer nPressureSaturation=size(ancillaryCoefficients.pressureSaturation, 1);
  Real[nPressureSaturation] n=ancillaryCoefficients.pressureSaturation[:, 1];
  Real[nPressureSaturation] theta=ancillaryCoefficients.pressureSaturation[:, 2];

  constant Real eps=1e-9;

algorithm
  assert(T <= T_crit+eps, "Ancillary.saturationPressure_T error: Temperature is higher than critical temperature", AssertionLevel.warning);
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
      <dd> <b>Eine mathematisch statistische Methode zum Aufstellen thermodynamischer Gleichungen - gezeigt am Beispiel der Dampfdruckkurve reiner fluider Stoffe.</b><br />
           Forschrittberichte der VDI Zeitschriften, Reihe 3, Nr. 39 (1974)
      </dd>
      </dl>
      </html>"));
end saturationPressure_T;
