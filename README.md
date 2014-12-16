# HelmholtzMedia
Modelica library for the calculation of fluid properties from an equation of state (EoS).

## Library description
This library calculates fluid properties from an equation of state (EoS) directly within Modelica. 
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.  

A general description of the library can be found in the related publications:
* Thorade, M. and Saadat, A. (2012). "[HelmholtzMedia - A fluid properties library][1]", 
_Proceedings of the 9th International Modelica Conference_, 
doi:10.3384/ecp1207663
* Thorade, M. (2012) "[Poster HelmholtzMedia][2]", 
_VDI Thermodynamik-Kolloquium 2012_  
* Thorade, M. and Saadat, A. (2013). "[Partial derivatives of thermodynamic state properties for dynamic simulation][3]", 
_Environmental Earth Sciences_, 
doi:10.1007/s12665-013-2394-z
* Thorade, M. (2014/unpublished). "[Entropiebasierte Bewertungskriterien für den Wärmeübergang in Kraftwerksprozessen und ihre Relevanz für praktische Anwendungen][4]", 
_Dissertation (TU Hamburg-Harburg)_

The following fluids have been implemented with EoS and transport properties:  
[Butane](HelmholtzMedia/HelmholtzFluids/Butane/package.mo), 
[Ethanol](HelmholtzMedia/HelmholtzFluids/Ethanol/package.mo), 
[Isobutane](HelmholtzMedia/HelmholtzFluids/Isobutane/package.mo), 
[Isopentane](HelmholtzMedia/HelmholtzFluids/Isopentane/package.mo), 
[Pentane](HelmholtzMedia/HelmholtzFluids/Pentane/package.mo), 
[Propane](HelmholtzMedia/HelmholtzFluids/Propane/package.mo), 
[R134a](HelmholtzMedia/HelmholtzFluids/R134a/package.mo) (with three reference states). 

The following fluids have been implemented with EoS only:  
[Helium](HelmholtzMedia/HelmholtzFluids/Helium/package.mo),
[Hexamethyldisiloxane (HMDS)](HelmholtzMedia/HelmholtzFluids/HMDS/package.mo). 

## Current release
Download the newest [tagged version](../../tags).  
In the future, there might be a release branch and official releases.

## License
Copyright © 2009-2013 Helmholtz Centre Potsdam, GFZ German Research Centre for Geosciences

This Modelica package is free software and the use is completely at your own risk;  
it can be redistributed and/or modified under the terms of the [Modelica License 2](https://www.modelica.org/licenses/ModelicaLicense2) or newer.  
For license conditions (including the disclaimer of warranty) visit [https://www.modelica.org/licenses/](https://www.modelica.org/licenses/).

## Development and contribution 
You may report feedback, issues or feature-requests using the [Issues](../../issues) button.  
Code contributions are very welcome, especially in the form of [Pull Requests](../../pulls).


[1]: http://goo.gl/7obgqq "Modelica Conference Paper: HelmholtzMedia implementation"
[2]: http://goo.gl/6dNHNP "Conference Poster: HelmholtzMedia implementation"
[3]: http://goo.gl/xbGJA9 "ISI Journal Paper: Partial derivatives"
[4]: http://goo.gl/4xgQ5X "Dissertation"
