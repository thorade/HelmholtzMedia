## HelmholtzFluids
This library is written in Modelica.
The purpose of this library is to calculate fluid properties from an equation of state (EoS), directly within Modelica and not from an external dll.
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.
I am planning to submit this library to the Modelica 2012 conference.

### Status
So far, the fluids n-Butane, R134a, Isobutane and Isopentane are implemented (but not fully validated)

### ToDo
* Remove optional inputs (make them protected variables instead) to get rid of warnings
* Improve BaseProperties (use preferredState or something else?)
* Add annotation regarding derivative and inverse functions
* Add annotation regarding Documentation (copy+paste+translate from conference paper)
* The next fluids to be implemented are probably 
  * propane, 
  * ammonia, 
  * CO2,
  * ethanol,
  * working fluids for Organic-Rankine-Cycles,
  * water.
* Low priority:
  * Add setState_Ts (needed for isotherms in h,s-diagram)
  * Add setState_Th (needed for numerical validation of (dd/dT)@h=const)