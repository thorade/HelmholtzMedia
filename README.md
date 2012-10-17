## HelmholtzMedia
This library is written in Modelica.
The purpose of this library is to calculate fluid properties from an equation of state (EoS), directly within Modelica and not from an external dll.
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.
Also see this [overview over the implementation][1].


### Implemented Fluids
* Butane
* Ethanol
* Isobutane
* Isopentane
* Pentane (hyperbolic terms by Martin Ryhl Kærn)
* Propane
* R134a

### Fluids that might be included in future versions
* Butane (short technical EoS)
* HFC-125 (should be easy)
* HFC-1234yf (ECS viscosity model not yet implemented)
* CO2 (non-analytical terms not yet implemented)
* water according to IAPWS-95 (non-analytical terms not yet implemented)
* nitrogen
* methane
* ammonia
* helium
* working fluids for Organic-Rankine-Cycles.

### Known issues & ToDo
* Index reduction doesn't work (there are some numerical Jacobians). Does this have to be fixed in `BaseProperties`, or in `density_ph`?
* Non-analytic critical terms for residual Helmholtz energy not yet implemented (needed for water and CO2)
* Reference state is fixed. Not sure how important choosing a differend reference state is.
* Add `setState_hs` (needed for turbine calculation, when power is given and p_out is to be determined)
* More testing would be nice (beta users welcome)
* Documentation could be extended

[1]: http://goo.gl/HeUzM "HelmholtzMedia CheatSheet"