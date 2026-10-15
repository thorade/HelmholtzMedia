within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Transport;
function dynamicViscosity_der "time derivative of dynamicViscosity"
  input Density d "density";
  input Temperature T "temperature";
  input FixedPhase phase=0 "phase of the state (not used)";
  input Real d_der "time derivative of density";
  input Real T_der "time derivative of temperature";
  output Real eta_der "time derivative of dynamic viscosity";

protected
  DynamicViscosityModel dynamicViscosityModel=dynamicViscosityCoefficients.dynamicViscosityModel;
  CollisionIntegralModel collisionIntegralModel=dynamicViscosityCoefficients.collisionIntegralModel;
  constant Real micro=1e-6;

  MolarMass MM=fluidConstants[1].molarMass;
  Real dm=d/(1000*MM) "molar density in mol/l";
  Real dm_d=1/(1000*MM) "derivative of dm w.r.t. d";
  Real epsilon_kappa=dynamicViscosityCoefficients.epsilon_kappa;
  Real sigma=dynamicViscosityCoefficients.sigma;

  // dilute contribution
  Real[size(dynamicViscosityCoefficients.a, 1),2] a=dynamicViscosityCoefficients.a;
  Real[size(dynamicViscosityCoefficients.CET, 1),2] CET=dynamicViscosityCoefficients.CET;
  Real[size(dynamicViscosityCoefficients.CET_den, 1),2] CET_den=dynamicViscosityCoefficients.CET_den;
  Temperature T_red_0=dynamicViscosityCoefficients.reducingTemperature_0;
  Real eta_red_0=dynamicViscosityCoefficients.reducingViscosity_0;
  Real T_star "reduced temperature of the collision integral";
  Real T_star_T "derivative of T_star w.r.t. T";
  Real S "sum in the collision integral";
  Real S_T "derivative of S w.r.t. T";
  Real Omega "collision integral";
  Real Omega_T "derivative of Omega w.r.t. T";
  Real tau "reduced temperature";
  Real num;
  Real num_tau;
  Real den;
  Real den_tau;
  Real eta_0 "dilute contribution (depends on T only)";
  Real eta_0_T;

  // initial density contribution
  Real[size(dynamicViscosityCoefficients.b, 1),2] b=dynamicViscosityCoefficients.b;
  Real eta_red_1=dynamicViscosityCoefficients.reducingViscosity_1;
  Real B "second viscosity virial coefficient, l/mol";
  Real B_T;
  Real eta_1_T;
  Real eta_1_d;

  // residual contribution
  Real[size(dynamicViscosityCoefficients.c, 1),1] c=dynamicViscosityCoefficients.c;
  Real[size(dynamicViscosityCoefficients.g, 1),2] g=dynamicViscosityCoefficients.g;
  Real[size(dynamicViscosityCoefficients.e, 1),5] e=dynamicViscosityCoefficients.e;
  Real[size(dynamicViscosityCoefficients.nu_po, 1),5] nu_po=dynamicViscosityCoefficients.nu_po;
  Real[size(dynamicViscosityCoefficients.de_po, 1),5] de_po=dynamicViscosityCoefficients.de_po;
  Temperature T_red_residual=dynamicViscosityCoefficients.reducingTemperature_residual;
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_residual=MM/dynamicViscosityCoefficients.reducingMolarVolume_residual;
  Real eta_red_residual=dynamicViscosityCoefficients.reducingViscosity_residual;
  Real tau_T "derivative of tau w.r.t. T";
  Real delta "reduced density";
  Real delta_d "derivative of delta w.r.t. d";
  Real delta_exp "reduced density in the exponential terms";
  Real delta_exp_d "derivative of delta_exp w.r.t. d";
  Real delta_0 "close packed density";
  Real delta_0_tau "derivative of delta_0 w.r.t. tau";
  Real p "polynomial term";
  Real p_T;
  Real p_d;
  Real X "exponential factor";
  Real X_d;
  Real xnum=0;
  Real xnum_T=0;
  Real xnum_d=0;
  Real xden=0;
  Real xden_T=0;
  Real xden_d=0;
  Real G;
  Real G_T;
  Real H;
  Real H_dm;
  Real F;
  Real F_T;
  Real F_dm;
  Real x "helium: log(T)";
  Real x_T;
  Real Bh;
  Real Bh_x;
  Real Ch;
  Real Ch_x;
  Real Dh;
  Real Dh_x;
  Real d_gcm=d/1000 "helium: density in g/cm3";
  Real E;
  Real eta_e;
  Real eta_e_T;
  Real eta_e_d;
  Real eta_0a;
  Real eta_0a_T;
  Real eta_0b;
  Real eta_0b_T;
  Real eta_0h;
  Real eta_0h_T;
  Real eta_r_T=0;
  Real eta_r_d=0;

algorithm
  // ---------------------------------------------------------------------------
  // residual contribution (dynamicViscosity_residual), partial derivatives w.r.t. T and d
  // ---------------------------------------------------------------------------
  if (dynamicViscosityModel == DynamicViscosityModel.VS0) then
    if mediumName == "helium" then
      if (T <= 300.0E0) then
        x := log(T);
        x_T := 1/T;
      else
        x := 5.703782474656201;
        x_T := 0;
      end if;
      Bh := -47.5295259E0/x + 87.6799309E0 - 42.0741589E0*x + 8.33128289E0*x^2 - 0.589252385E0*x^3;
      Ch :=  547.309267E0/x - 904.870586E0 + 431.404928E0*x - 81.4504854E0*x^2 + 5.37008433E0*x^3;
      Dh := -1684.39324E0/x + 3331.08630E0 - 1632.19172E0*x + 308.804413E0*x^2 - 20.2936367E0*x^3;
      Bh_x := 47.5295259E0/x^2 - 42.0741589E0 + 2*8.33128289E0*x - 3*0.589252385E0*x^2;
      Ch_x := -547.309267E0/x^2 + 431.404928E0 - 2*81.4504854E0*x + 3*5.37008433E0*x^2;
      Dh_x := 1684.39324E0/x^2 - 1632.19172E0 + 2*308.804413E0*x - 3*20.2936367E0*x^2;
      E := Bh*d_gcm + Ch*d_gcm^2 + Dh*d_gcm^3;
      eta_e := exp(min(100, E));
      if (E < 100) then
        eta_e_T := eta_e*(Bh_x*d_gcm + Ch_x*d_gcm^2 + Dh_x*d_gcm^3)*x_T;
        eta_e_d := eta_e*(Bh + 2*Ch*d_gcm + 3*Dh*d_gcm^2)/1000;
      else
        eta_e_T := 0;
        eta_e_d := 0;
      end if;
      eta_0a := exp(-0.135311743E0/x + 1.00347841E0 + 1.20654649E0*x - 0.149564551E0*x^2 + 0.0125208416E0*x^3);
      eta_0a_T := eta_0a*(0.135311743E0/x^2 + 1.20654649E0 - 2*0.149564551E0*x + 3*0.0125208416E0*x^2)*x_T;
      if (T < 100.0E0) then
        // eta_r = eta_0a*eta_e
        eta_r_T := eta_0a_T*eta_e + eta_0a*eta_e_T;
        eta_r_d := eta_0a*eta_e_d;
      else
        eta_0b := 196.0E0*T^0.71938E0*exp(12.451E0/T - 295.67E0/T^2 - 4.1249E0);
        eta_0b_T := eta_0b*(0.71938E0/T - 12.451E0/T^2 + 2*295.67E0/T^3);
        if (T <= 110.0E0) then
          eta_0h := eta_0a + (eta_0b - eta_0a)*(T - 100.0E0)/10.0E0;
          eta_0h_T := eta_0a_T + (eta_0b_T - eta_0a_T)*(T - 100.0E0)/10.0E0 + (eta_0b - eta_0a)/10.0E0;
        else
          eta_0h := eta_0b;
          eta_0h_T := eta_0b_T;
        end if;
        // eta_r = eta_0a*eta_e + eta_0 - eta_0a
        eta_r_T := eta_0a_T*eta_e + eta_0a*eta_e_T + eta_0h_T - eta_0a_T;
        eta_r_d := eta_0a*eta_e_d;
      end if;
      eta_r_T := eta_r_T/10.0E0;
      eta_r_d := eta_r_d/10.0E0;
    else
      assert(false, "dynamicViscosity_der: no hardcoded VS0 model for " + mediumName);
    end if;

  elseif ((dynamicViscosityModel == DynamicViscosityModel.VS1)
      or  (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)) then
    tau := T/T_red_residual;
    tau_T := 1/T_red_residual;
    delta := d/d_red_residual;
    delta_d := 1/d_red_residual;
    if (abs(d_red_residual - 1) > 0.001) then
      delta_exp := d/d_crit;
      delta_exp_d := 1/d_crit;
    else
      delta_exp := delta;
      delta_exp_d := delta_d;
    end if;

    if size(g, 1) > 0 then
      if (dynamicViscosityModel == DynamicViscosityModel.VS1) then
        delta_0 := sum(g[i, 1]*tau^g[i, 2] for i in 1:size(g, 1));
        delta_0_tau := sum(if g[i, 2] <> 0 then g[i, 1]*g[i, 2]*tau^(g[i, 2] - 1) else 0 for i in 1:size(g, 1));
      else
        // VS1_alternative: delta_0 = g1/(1 + S)
        S := sum(g[i, 1]*tau^g[i, 2] for i in 2:size(g, 1));
        delta_0 := g[1, 1]/(1 + S);
        delta_0_tau := -g[1, 1]/(1 + S)^2*sum(if g[i, 2] <> 0 then g[i, 1]*g[i, 2]*tau^(g[i, 2] - 1) else 0 for i in 2:size(g, 1));
      end if;
    else
      delta_0 := 1.0;
      delta_0_tau := 0;
    end if;

    // simple polynomial terms, optionally with exponential factor
    for i in 1:size(e, 1) loop
      p := e[i, 1]*tau^e[i, 2]*delta^e[i, 3]*delta_0^e[i, 4];
      p_T := e[i, 1]*delta^e[i, 3]*((if e[i, 2] <> 0 then e[i, 2]*tau^(e[i, 2] - 1)*delta_0^e[i, 4] else 0)
           + (if e[i, 4] <> 0 then tau^e[i, 2]*e[i, 4]*delta_0^(e[i, 4] - 1)*delta_0_tau else 0))*tau_T;
      p_d := if e[i, 3] <> 0 then e[i, 1]*tau^e[i, 2]*e[i, 3]*delta^(e[i, 3] - 1)*delta_0^e[i, 4]*delta_d else 0;
      if (e[i, 5] > 0) then
        X := exp(-delta_exp^e[i, 5]);
        X_d := -X*e[i, 5]*delta_exp^(e[i, 5] - 1)*delta_exp_d;
        eta_r_T := eta_r_T + p_T*X;
        eta_r_d := eta_r_d + p_d*X + p*X_d;
      else
        eta_r_T := eta_r_T + p_T;
        eta_r_d := eta_r_d + p_d;
      end if;
    end for;

    // rational polynomial terms (RefProp algorithm: the running sum is multiplied by the exponential factor)
    for i in 1:size(nu_po, 1) loop
      p := nu_po[i, 1]*tau^nu_po[i, 2]*delta^nu_po[i, 3]*delta_0^nu_po[i, 4];
      p_T := nu_po[i, 1]*delta^nu_po[i, 3]*((if nu_po[i, 2] <> 0 then nu_po[i, 2]*tau^(nu_po[i, 2] - 1)*delta_0^nu_po[i, 4] else 0)
           + (if nu_po[i, 4] <> 0 then tau^nu_po[i, 2]*nu_po[i, 4]*delta_0^(nu_po[i, 4] - 1)*delta_0_tau else 0))*tau_T;
      p_d := if nu_po[i, 3] <> 0 then nu_po[i, 1]*tau^nu_po[i, 2]*nu_po[i, 3]*delta^(nu_po[i, 3] - 1)*delta_0^nu_po[i, 4]*delta_d else 0;
      xnum := xnum + p;
      xnum_T := xnum_T + p_T;
      xnum_d := xnum_d + p_d;
      if (nu_po[i, 5] > 0) then
        X := exp(-delta_exp^nu_po[i, 5]);
        X_d := -X*nu_po[i, 5]*delta_exp^(nu_po[i, 5] - 1)*delta_exp_d;
        xnum_d := xnum_d*X + xnum*X_d;
        xnum_T := xnum_T*X;
        xnum := xnum*X;
      end if;
    end for;
    for i in 1:size(de_po, 1) loop
      p := de_po[i, 1]*tau^de_po[i, 2]*delta^de_po[i, 3]*delta_0^de_po[i, 4];
      p_T := de_po[i, 1]*delta^de_po[i, 3]*((if de_po[i, 2] <> 0 then de_po[i, 2]*tau^(de_po[i, 2] - 1)*delta_0^de_po[i, 4] else 0)
           + (if de_po[i, 4] <> 0 then tau^de_po[i, 2]*de_po[i, 4]*delta_0^(de_po[i, 4] - 1)*delta_0_tau else 0))*tau_T;
      p_d := if de_po[i, 3] <> 0 then de_po[i, 1]*tau^de_po[i, 2]*de_po[i, 3]*delta^(de_po[i, 3] - 1)*delta_0^de_po[i, 4]*delta_d else 0;
      xden := xden + p;
      xden_T := xden_T + p_T;
      xden_d := xden_d + p_d;
      if (de_po[i, 5] > 0) then
        X := exp(-delta_exp^de_po[i, 5]);
        X_d := -X*de_po[i, 5]*delta_exp^(de_po[i, 5] - 1)*delta_exp_d;
        xden_d := xden_d*X + xden*X_d;
        xden_T := xden_T*X;
        xden := xden*X;
      end if;
    end for;
    if xden <> 0 then
      eta_r_T := eta_r_T + (xnum_T*xden - xnum*xden_T)/xden^2;
      eta_r_d := eta_r_d + (xnum_d*xden - xnum*xden_d)/xden^2;
    end if;

  elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
    G := c[1, 1] + c[2, 1]/T;
    G_T := -c[2, 1]/T^2;
    H := sqrt(dm)*(dm - c[8, 1])/c[8, 1];
    H_dm := (1.5*dm - 0.5*c[8, 1])/(c[8, 1]*sqrt(dm));
    F := G + (c[3, 1] + c[4, 1]*T^(-3/2))*dm^0.1 + (c[5, 1] + c[6, 1]/T + c[7, 1]/T^2)*H;
    F_T := G_T - 1.5*c[4, 1]*T^(-5/2)*dm^0.1 - (c[6, 1]/T^2 + 2*c[7, 1]/T^3)*H;
    F_dm := (c[3, 1] + c[4, 1]*T^(-3/2))*0.1*dm^(-0.9) + (c[5, 1] + c[6, 1]/T + c[7, 1]/T^2)*H_dm;
    // eta_r = exp(F) - exp(G)
    eta_r_T := exp(F)*F_T - exp(G)*G_T;
    eta_r_d := exp(F)*F_dm*dm_d;
  else
    assert(false, "dynamicViscosity_der: unknown or not implemented dynamicViscosityModel");
  end if;
  eta_r_T := eta_r_T*eta_red_residual;
  eta_r_d := eta_r_d*eta_red_residual;

  if (dynamicViscosityModel == DynamicViscosityModel.VS0) then
    // hardcoded models: residual contribution is the full viscosity
    eta_der := micro*(eta_r_T*T_der + eta_r_d*d_der);
  else
    // -------------------------------------------------------------------------
    // collision integral and dilute contribution (dynamicViscosity_dilute), depends on T only
    // -------------------------------------------------------------------------
    if (collisionIntegralModel == CollisionIntegralModel.CI0) then
      T_star := T/epsilon_kappa;
      T_star_T := 1/epsilon_kappa;
      Omega := 1.16145/T_star^0.14874 + 0.52487*exp(-0.77320*T_star) + 2.16178*exp(-2.43787*T_star);
      Omega_T := (-0.14874*1.16145/T_star^1.14874 - 0.77320*0.52487*exp(-0.77320*T_star)
                 - 2.43787*2.16178*exp(-2.43787*T_star))*T_star_T;
    elseif (collisionIntegralModel == CollisionIntegralModel.CI1) then
      T_star := Modelica.Math.log(T/epsilon_kappa);
      if abs(T_star) < 1E-20 then
        T_star := 1E-20;
      end if;
      T_star_T := 1/T;
      Omega := exp(sum(a[i, 1]*T_star^a[i, 2] for i in 1:size(a, 1)));
      Omega_T := Omega*sum(if a[i, 2] <> 0 then a[i, 1]*a[i, 2]*T_star^(a[i, 2] - 1) else 0 for i in 1:size(a, 1))*T_star_T;
    elseif (collisionIntegralModel == CollisionIntegralModel.CI2) then
      T_star := (epsilon_kappa/T)^(1/3);
      T_star_T := -T_star/(3*T);
      S := sum(a[i, 1]*T_star^(4 - i) for i in 1:size(a, 1));
      S_T := sum(a[i, 1]*(4 - i)*T_star^(3 - i) for i in 1:size(a, 1))*T_star_T;
      Omega := 1/S;
      Omega_T := -S_T/S^2;
    else
      assert(false, "dynamicViscosity_der: unknown CollisionIntegralModel");
    end if;

    if ((dynamicViscosityModel == DynamicViscosityModel.VS1)
     or (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)) then
      tau := T/T_red_0;
      if size(CET_den, 1) > 0 then
        num := sum(CET[i, 1]*tau^CET[i, 2] for i in 1:size(CET, 1));
        num_tau := sum(if CET[i, 2] <> 0 then CET[i, 1]*CET[i, 2]*tau^(CET[i, 2] - 1) else 0 for i in 1:size(CET, 1));
        den := sum(CET_den[i, 1]*tau^CET_den[i, 2] for i in 1:size(CET_den, 1));
        den_tau := sum(if CET_den[i, 2] <> 0 then CET_den[i, 1]*CET_den[i, 2]*tau^(CET_den[i, 2] - 1) else 0 for i in 1:size(CET_den, 1));
        eta_0 := num/den;
        eta_0_T := (num_tau*den - num*den_tau)/den^2/T_red_0;
      else
        eta_0 := CET[1, 1]*sqrt(tau)/(sigma^2*Omega)
               + sum(CET[i, 1]*tau^CET[i, 2] for i in 2:size(CET, 1));
        eta_0_T := CET[1, 1]/sigma^2*(0.5/(sqrt(tau)*T_red_0*Omega) - sqrt(tau)*Omega_T/Omega^2)
                 + sum(if CET[i, 2] <> 0 then CET[i, 1]*CET[i, 2]*tau^(CET[i, 2] - 1) else 0 for i in 2:size(CET, 1))/T_red_0;
      end if;
    elseif (dynamicViscosityModel == DynamicViscosityModel.VS2) then
      eta_0 := CET[1, 1]*T^CET[1, 2]/(sigma^2*Omega);
      eta_0_T := eta_0*(CET[1, 2]/T - Omega_T/Omega);
    else
      assert(false, "dynamicViscosity_der: unknown or not implemented dynamicViscosityModel");
    end if;
    eta_0 := eta_0*eta_red_0;
    eta_0_T := eta_0_T*eta_red_0;

    // -------------------------------------------------------------------------
    // initial density contribution (dynamicViscosity_initial)
    // -------------------------------------------------------------------------
    if ((dynamicViscosityModel == DynamicViscosityModel.VS1)
     or (dynamicViscosityModel == DynamicViscosityModel.VS1_alternative)) then
      // eta_1 = eta_0*B*dm (Rainwater and Friend)
      T_star := T/epsilon_kappa;
      B := sum(b[i, 1]*T_star^b[i, 2] for i in 1:size(b, 1))*0.6022137*sigma^3;
      B_T := sum(if b[i, 2] <> 0 then b[i, 1]*b[i, 2]*T_star^(b[i, 2] - 1) else 0 for i in 1:size(b, 1))*0.6022137*sigma^3/epsilon_kappa;
      eta_1_T := (eta_0_T*B + eta_0*B_T)*dm;
      eta_1_d := eta_0*B*dm_d;
    else
      // VS2: eta_1 = dm*(b1 + b2*(b3 - log(T/b4))^2)
      eta_1_T := -dm*b[2, 1]*2*(b[3, 1] - log(T/b[4, 1]))/T;
      eta_1_d := (b[1, 1] + b[2, 1]*(b[3, 1] - log(T/b[4, 1]))^2)*dm_d;
    end if;
    eta_1_T := eta_1_T*eta_red_1;
    eta_1_d := eta_1_d*eta_red_1;

    eta_der := micro*((eta_0_T + eta_1_T + eta_r_T)*T_der + (eta_1_d + eta_r_d)*d_der);
  end if;

  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Transport.dynamicViscosity\">dynamicViscosity</a>:
the partial derivatives of the dilute, initial density and residual contributions with respect to T and d are
computed analytically in one function (the dilute contribution and the collision integral are evaluated once and
shared with the initial density contribution), then
eta_der = deta/dT*T_der + deta/dd*d_der.
The model branches are the same as in dynamicViscosity_dilute, dynamicViscosity_initial and dynamicViscosity_residual.
</p>
</html>"));
end dynamicViscosity_der;
