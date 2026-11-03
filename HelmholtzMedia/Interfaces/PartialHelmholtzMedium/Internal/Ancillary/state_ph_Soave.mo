within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Ancillary;
function state_ph_Soave
  "ancillary function: approximate temperature and density of vapour and gas-like states for given p and h (ideal gas plus Redlich-Kwong-Soave departure), used as start values"

  input Modelica.Units.SI.AbsolutePressure p;
  input Modelica.Units.SI.SpecificEnthalpy h;
  input Modelica.Units.SI.Temperature T_low "lower bound of T, e.g. the saturation temperature";
  input Modelica.Units.SI.Temperature T_start "start value of the iteration";
  output Modelica.Units.SI.Temperature T;
  output Modelica.Units.SI.Density d;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM "specific gas constant";
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
  constant Temperature T_reducing=helmholtzCoefficients.T_reducing "reducing temperature of the EoS";
  constant Temperature T_high=1.5*fluidLimits.TMAX;

  // RKS: Redlich-Kwong-Soave, same constants as in density_pT_Soave
  constant Real omega = fluidConstants[1].acentricFactor;
  constant Real m = 0.480 + 1.574*omega - 0.176*omega^2;
  constant Real Omega_a = 1/(9*(2^(1/3) - 1));
  constant Real Omega_b = (2^(1/3) - 1)/3;
  constant Real a_crit = Omega_a*R_s^2*T_crit^2/p_crit;
  constant Real b = Omega_b*R_s*T_crit/p_crit;

  constant Real tolerance=1e-4 "relative tolerance of T (start value only)";
  constant Integer iter_max=20;

  Real tau;
  Real f_i;
  Real f_it;
  Real f_itt;
  Real f_ittt;
  SpecificEnthalpy h0 "ideal-gas enthalpy";
  SpecificHeatCapacity cp0 "ideal-gas heat capacity";
  Real sa "square root of the Soave alpha function";
  Real a;
  Real dadT;
  Real Z;
  SpecificEnthalpy h_dep "departure function h - h0 of RKS";
  SpecificEnthalpy h_dep_old = 0;
  Temperature T_old = 0;
  Real slope;
  Temperature dT = 1;
  Integer iter = 0;

algorithm
  T := min(max(T_start, T_low), T_high);
  d := 0;
  while (iter < iter_max) and (abs(dT) > tolerance*T) loop
    iter := iter + 1;
    // ideal part of the EoS itself, so that h0 includes the reference state of the medium
    tau := T_reducing/T;
    (f_i, f_it, f_itt, f_ittt) := EoS.f_i_derivs(delta=1, tau=tau, order=2);
    h0 := R_s*T*(1 + tau*f_it);
    cp0 := R_s*(1 - tau^2*f_itt);
    // departure from RKS, vapour root
    d := density_pT_Soave(T=T, p=p, psat=2*p);
    Z := p/(d*R_s*T);
    sa := 1 + m*(1 - sqrt(T/T_crit));
    a := a_crit*sa^2;
    dadT := -a_crit*m*sa/sqrt(T*T_crit);
    h_dep := R_s*T*(Z - 1) + (T*dadT - a)/b*log(1 + b*d);
    // Newton step on h0(T) + h_dep(T) = h; the temperature derivative of h_dep is taken from the previous step
    slope := cp0 + (if (iter > 1) and (abs(T - T_old) > 0) then (h_dep - h_dep_old)/(T - T_old) else 0);
    slope := max(slope, 0.5*cp0);
    T_old := T;
    h_dep_old := h_dep;
    dT := (h - h0 - h_dep)/slope;
    T := min(max(T + dT, T_low), T_high);
  end while;
  d := density_pT_Soave(T=T, p=p, psat=2*p);

  annotation (Documentation(info="<html>
<p>
Start values for vapour and gas-like states in <code>setState_phX</code>. The temperature solves
h<sup>0</sup>(T) + h<sup>dep</sup>(T, p) = h, where h<sup>0</sup> is the ideal-gas enthalpy from the ideal part of the
Helmholtz EoS (including the reference state) and h<sup>dep</sup> the departure function of the Redlich-Kwong-Soave
equation (Soave 1972); the density is the vapour root of the Redlich-Kwong-Soave equation at this temperature.
The function needs no evaluation of the residual part of the EoS.
</p>
</html>"));
end state_ph_Soave;
