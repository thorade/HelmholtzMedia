within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function density_pT_Soave "solve Redlich-Kwong-Soave for density"
  input Temperature T;
  input AbsolutePressure p;
  input AbsolutePressure psat=Ancillary.saturationPressure_T(T=T);
  output Density d;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // RKS: Redlich-Kwong-Soave (see Soave 1979)
  constant Real omega = fluidConstants[1].acentricFactor;
  constant Real m = 0.480 + 1.574*omega - 0.176*omega^2;
  Real a = 0.42747*R^2*T_crit^2/p_crit*(1 + m*(1 - sqrt(T/T_crit)))^2;
  Real b = 0.08664*R*T_crit/p_crit;
  Real A = a*p/(R^2*T^2);
  Real B = b*p/(R*T);
  Real r = (A-B-B^2)-1/3;
  Real q = -2/27 + 1/3*(A-B-B^2) - A*B;
  Real D = (r/3)^3 + (q/2)^2 "discriminant";
  Real u;
  Real u3;
  Real Y1;
  Real Y2;
  Real Y3;
  Real Theta;
  Real phi;
  import Modelica.Constants.pi;

algorithm
  // get density start value from Redlich-Kwong-Soave (see Span 2000, section "3.3.1 Calculations based on pT" )
  // Modelica.Utilities.Streams.print("", "printlog.txt");
  // Modelica.Utilities.Streams.print("Redlich-Kwong-Soave called with p=" + String(p) + " and T=" +String(T), "printlog.txt");
  // Modelica.Utilities.Streams.print("RKS discriminant D=" + String(D) + " and r=" + String(r), "printlog.txt");
  if (D >= 0) then
    u3 := -q/2 + sqrt(D);
    u := sign(u3)*abs(u3)^(1/3);
    Y1 := u-r/(3*u);
    // Modelica.Utilities.Streams.print("RKS has one root (Y1=" + String(Y1) + ")", "printlog.txt");
    d := p/(R*T*(Y1+1/3));
  elseif ((abs(D) < 1e-8) and (abs(r) < 1e-3)) then
    // Modelica.Utilities.Streams.print("close to critical region, use critical density");
    d := d_crit;
  else
    // case D<0
    Theta := sqrt(-r^3/27);
    phi := acos(-q/(2*Theta));
    Y1 := 2*Theta^(1/3)*cos(phi/3);
    Y2 := 2*Theta^(1/3)*cos(phi/3+2*pi/3);
    Y3 := 2*Theta^(1/3)*cos(phi/3+4*pi/3);
    // Modelica.Utilities.Streams.print("RKS has three possible roots(Y1=" + String(Y1) + ", Y2=" + String(Y2) + ", Y3=" + String(Y3), "printlog.txt");
    if (T <= T_crit) then
      // Modelica.Utilities.Streams.print("T<T_crit: multiple roots due to phase boundary", "printlog.txt");
      // Modelica.Utilities.Streams.print("d(Y1)=" + String(p/(R*T*(Y1+1/3))) + ", d(Y2)=" + String(p/(R*T*(Y2+1/3))) + ", d(p/(R*T*(Y3+1/3)))=" + String(Y3), "printlog.txt");
      if (p > psat) then
        // liquid, use smallest Y
        Y1 := min(Y1,Y2);
        Y1 := min(Y1,Y3);
        d := p/(R*T*(Y1+1/3));
      elseif (p < psat) then
        // vapor, use largest Y
        Y1 := max(Y1,Y2);
        Y1 := max(Y1,Y3);
        d := p/(R*T*(Y1+1/3));
      else
        // this is very unlikely, but not impossible
        assert(p <> psat, "setState_pTX error: pressure equals saturation pressure");
      end if;
    else
      // Modelica.Utilities.Streams.print("T>T_crit: multiple roots can occur, but two of the roots result in negative densities", "printlog.txt");
      // Modelica.Utilities.Streams.print("d(Y1)=" + String(p/(R*T*(Y1+1/3))) + ", d(Y2)=" + String(p/(R*T*(Y2+1/3))) + ", d(Y3)=" + String(p/(R*T*(Y3+1/3))), "printlog.txt");
      d := max(p/(R*T*(Y1+1/3)), p/(R*T*(Y2+1/3)));
      d := max(p/(R*T*(Y3+1/3)), d);
    end if;
  end if;
  // Modelica.Utilities.Streams.print("RKS finished, d=" + String(d), "printlog.txt");

end density_pT_Soave;
