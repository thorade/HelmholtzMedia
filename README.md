## HelmholtzFluids
This library is written in Modelica.
The purpose of this library is to calculate fluid properties from an equation of state (EoS), directly within Modelica and not from an external dll.
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.
I am planning to submit this library to the Modelica 2012 conference.

### Status
So far, the fluids n-Butane, R134a, Isobutane, Isopentane and Ethanol are implemented (but not fully validated)

### ToDo
* Replace asserts in setSat_T with extrapolation
* Improve `BaseProperties` (should independentVariables be used here?)
* Check which functions from `TwoPhaseMedium` are still missing & implement them
* Add `annotation(inverse(b=b(a))` where possible & not yet implemented
* Add `annotation(derivative=...)` to `temperature_ph` ~~and `density_ph`~~
* Add `annotation(Documentation(info=...html...))` using copy+paste from conference paper
* The next fluids to be implemented are probably 
  * propane, 
  * ammonia,
  * working fluids for Organic-Rankine-Cycles.
* Add non-analytic critical terms for residual Helmholtz energy (needed for water and CO2)
* Add hyperbolic terms for ideal Helmholtz energy (used in short technical EoS and GERG-2008)
* Add choice of reference state
* Low priority:
  * Add `setState_Ts` (needed for isotherms in h,s-diagram)
  * Add `setState_Th` (needed for numerical validation of (dd/dT)@h=const)