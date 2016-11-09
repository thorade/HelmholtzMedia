% manually delete ReferenceState.csv
% uncomment print header in ReferenceState.mo
% set idealPower[1,1] and idealPower[2,1] to 0
% run ReferenceState.mo
% comment print header
% change idealPower[1,1] and [2,1] to some value different from 0
% run ReferenceState.mo several times using different values
% ------------------------------------------------------------------------
% start Matlab
% open ReferenceState.csv
% import as column vectors
% run ReferenceState.m
% ------------------------------------------------------------------------
% copy root1 to idealPower[1,1] and root2 to idealPower[2,1]
% ------------------------------------------------------------------------
% y=m*x+b
% for y=0 x=-b/m
% for y=200
% ------------------------------------------------------------------------

% Entropy
% find linear fit for reverse function
[polyS,StructureS,muS] = polyfit(sref,idealPower1,1)
% evaluate reverse function at desired value
sref0=1e3;
[idealPower1_ref,deltaS] = polyval(polyS,sref0,StructureS,muS)


% Enthalpy
% find linear fit for reverse function
[polyH,StructureH,muH] = polyfit(href,idealPower2,1)
% evaluate reverse function at desired href value
href0=2e5;
[idealPower2_ref,deltaH] = polyval(polyH,href0,StructureH,muH)
