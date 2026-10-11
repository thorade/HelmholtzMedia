within HelmholtzMedia.Examples.ConvergenceTest;
model Ancillary_density_pT_Soave
  "check the Redlich-Kwong-Soave start values for density against the EoS in the whole validity range"
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  parameter Integer nT = 20 "number of temperatures between TMIN and TMAX";
  parameter Integer np = 20 "number of pressures between p_low and PMAX, logarithmically spaced";
  parameter Real tol_vap = 0.1 "tolerance for abs(log(d_RKS/d_EoS)) in the vapour region";
  parameter Real tol_dense = 1.0 "tolerance for abs(log(d_RKS/d_EoS)) in the liquid and supercritical region";
  Real maxErr_vap(start=0, fixed=true) "largest abs(log(d_RKS/d_EoS)), vapour";
  Real maxErr_liq(start=0, fixed=true) "largest abs(log(d_RKS/d_EoS)), liquid";
  Real maxErr_sup(start=0, fixed=true) "largest abs(log(d_RKS/d_EoS)), supercritical";
  Integer nFailed(start=0, fixed=true) "number of failed checks";

protected
  constant Medium.Temperature T_crit = Medium.fluidConstants[1].criticalTemperature;
  constant Medium.AbsolutePressure p_crit = Medium.fluidConstants[1].criticalPressure;
  constant Medium.Temperature TMIN = Medium.fluidLimits.TMIN;
  constant Medium.Temperature TMAX = Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure PMAX = Medium.fluidLimits.PMAX;
  constant Medium.AbsolutePressure p_low = 1e-6;
  constant Medium.SpecificHeatCapacity R_s = Modelica.Constants.R/Medium.fluidConstants[1].molarMass;
  constant Medium.Density d_cov = p_crit/(0.08664*R_s*T_crit) "inverse RKS covolume 1/b";
  final constant String fileName = "soave.txt";
  Medium.Temperature T;
  Medium.AbsolutePressure p;
  Medium.AbsolutePressure psat;
  Medium.Density d_EoS;
  Medium.Density d_RKS;
  Real err;
  Integer region "1 vapour, 2 liquid, 3 supercritical";

algorithm
  when initial() then
    Modelica.Utilities.Files.remove(fileName);
    for i in 1:nT loop
      T := TMIN + (TMAX - TMIN)*(i - 1)/(nT - 1);
      psat := if (T < T_crit) then Medium.saturationPressure(T) else p_crit;
      for j in 1:np loop
        p := p_low*(PMAX/p_low)^((j - 1)/(np - 1));
        region := if (T >= T_crit) then 3 elseif (p < psat) then 1 else 2;
        // skip states close to the saturation curve, where the phase of the RKS root is ambiguous
        if (region == 3) or (abs(p/psat - 1) > 0.02) then
          d_EoS := Medium.density(Medium.setState_pTX(p=p, T=T, phase=1));
          d_RKS := Medium.Ancillary.density_pT_Soave(T=T, p=p, psat=psat);
          err := if (d_RKS > 0) then abs(log(d_RKS/d_EoS)) else Modelica.Constants.inf;
          if (region == 1) then
            maxErr_vap := max(maxErr_vap, err);
          elseif (region == 2) then
            maxErr_liq := max(maxErr_liq, err);
          else
            maxErr_sup := max(maxErr_sup, err);
          end if;
          if (d_RKS <= 0) or (d_RKS >= d_cov) or (err > (if (region == 1) then tol_vap else tol_dense)) then
            nFailed := nFailed + 1;
            Modelica.Utilities.Streams.print("Ancillary_density_pT_Soave FAILED: T=" + String(T) + " p=" + String(p)
              + " region=" + String(region) + " d_EoS=" + String(d_EoS) + " d_RKS=" + String(d_RKS));
            Modelica.Utilities.Streams.print("FAILED: T=" + String(T) + " p=" + String(p) + " region=" + String(region)
              + " d_EoS=" + String(d_EoS) + " d_RKS=" + String(d_RKS), fileName);
          end if;
        end if;
      end for;
    end for;
    Modelica.Utilities.Streams.print("maxErr_vap=" + String(maxErr_vap) + " maxErr_liq=" + String(maxErr_liq)
      + " maxErr_sup=" + String(maxErr_sup) + " nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, String(nFailed) + " RKS start value(s) failed for " + Medium.fluidConstants[1].chemicalFormula
      + ", see the lines starting with Ancillary_density_pT_Soave FAILED in the simulation log, or " + fileName);
  end when;

  annotation (Documentation(info="<html>
<p>Compares the density from the Redlich-Kwong-Soave equation of state (<code>Ancillary.density_pT_Soave</code>),
which is used as start value in <code>setState_pTX</code>, with the density from the Helmholtz EoS
on a grid of nT temperatures between TMIN and TMAX and np pressures between 1e-6 Pa and PMAX (logarithmically spaced).
States within 2 % of the saturation pressure are skipped.</p>
<p>The RKS density must be positive, smaller than the inverse covolume 1/b, and close to the EoS density:
abs(log(d_RKS/d_EoS)) &le; tol_vap in the vapour region and &le; tol_dense in the liquid and supercritical region,
where the cubic equation is only qualitatively correct. The largest deviations per region are available as
variables and, together with all failed states, written to <code>soave.txt</code>.
Use <code>redeclare package Medium</code> to check other fluids.</p>
<p>The test points include the case from issue #3 (Pentane at 1e-6 Pa and high temperature),
where the previous implementation picked a root with a specific volume smaller than the covolume.</p>
</html>"));
end Ancillary_density_pT_Soave;
