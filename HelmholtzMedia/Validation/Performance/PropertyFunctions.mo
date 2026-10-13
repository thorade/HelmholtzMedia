within HelmholtzMedia.Validation.Performance;
model PropertyFunctions "Performance test: repeated calls of setState and property functions"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Propane;
  parameter Integer nT=20 "number of temperatures";
  parameter Integer np=20 "number of pressures";
  parameter Integer nRepeat=5 "number of repetitions of the whole grid";
  parameter Boolean transport=true "include viscosity and thermal conductivity (set false for fluids without transport properties)";
  parameter String fileName="perf_PropertyFunctions.txt";
  Real checksum(unit="1") "sum of all computed properties, each divided by a reference value";
  Integer nCalls "number of grid points evaluated";

protected
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.AbsolutePressure p_crit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.Density d_ref=1;
  constant Medium.SpecificHeatCapacity cp_ref=1000;
  constant Medium.VelocityOfSound c_ref=1000;
  constant Medium.DynamicViscosity eta_ref=1e-6;
  constant Medium.ThermalConductivity lambda_ref=1e-3;
  constant Medium.Temperature T_ref=1;
  constant Medium.Temperature T_min=max(Medium.fluidLimits.TMIN*1.05, 0.6*T_crit);
  constant Medium.Temperature T_max=min(Medium.fluidLimits.TMAX, 2.5*T_crit);
  constant Medium.AbsolutePressure p_min=1e-3*p_crit;
  constant Medium.AbsolutePressure p_max=min(Medium.fluidLimits.PMAX, 5*p_crit);
  Medium.Temperature T;
  Medium.AbsolutePressure p;
  Medium.AbsolutePressure psat;
  Medium.ThermodynamicState state;
  Medium.ThermodynamicState state_ph;
  Medium.SaturationProperties sat;

algorithm
  when initial() then
    checksum := 0;
    nCalls := 0;
    for k in 1:nRepeat loop
      for i in 1:nT loop
        T := T_min*(T_max/T_min)^((i - 1)/(nT - 1));
        if T < 0.99*T_crit then
          sat := Medium.setSat_T(T=T);
          psat := sat.psat;
          checksum := checksum + sat.psat/p_crit + sat.liq.d/d_ref + sat.vap.d/d_ref;
        else
          psat := -1;
        end if;
        for j in 1:np loop
          p := p_min*(p_max/p_min)^((j - 1)/(np - 1));
          if (psat < 0) or (abs(p/psat - 1) > 0.02) then
            state := Medium.setState_pTX(p=p, T=T, X=Medium.reference_X);
            state_ph := Medium.setState_phX(p=p, h=state.h, X=Medium.reference_X);
            checksum := checksum + state.d/d_ref + state_ph.T/T_ref + Medium.specificHeatCapacityCp(state)/cp_ref
              + Medium.velocityOfSound(state)/c_ref;
            if transport then
              checksum := checksum + Medium.dynamicViscosity(state)/eta_ref + Medium.thermalConductivity(state)/lambda_ref;
            end if;
            nCalls := nCalls + 1;
          end if;
        end for;
      end for;
    end for;
    Modelica.Utilities.Files.remove(fileName);
    Modelica.Utilities.Streams.print("nCalls=" + String(nCalls) + " checksum=" + String(checksum, significantDigits=17), fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>
On a grid of temperatures (0.6*T_crit to 2.5*T_crit) and pressures (1e-3*p_crit to 5*p_crit), this model calls
setSat_T, setState_pTX, setState_phX, cp, speed of sound and, optionally, viscosity and thermal conductivity, nRepeat times.
Points within 2&nbsp;% of the saturation pressure are skipped. The checksum allows to verify that a different implementation gives the same results.
</p>
</html>"));
end PropertyFunctions;
