% manually delete ReferenceState.csv
% uncomment print header in ReferenceState.mo
% set idealPower[1,1] and idealPower[2,1] to 0 
% run ReferenceState.mo
% comment print header
% change idealPower[1,1] and [2,1] to some value different from 0
% run ReferenceState.mo several times using different values
% start Matlab
% open ReferenceState.csv
% import as column vectors
% run ReferenceState.m
% copy root1 to idealPower[1,1] and root2 to idealPower[2,1]

% y=m*x+b
% for y=0 x=-b/m

b1 = sref(1, 1);
length1 = size(sref,1);
dy1 =        sref(length1,1) -        sref(1,1);
dx1 = idealPower1(length1,1) - idealPower1(1,1);
m1= dy1/dx1;
root1 = -b1/m1


b2 = href(1,1);
length2 = size(href,1);
dy2 =        href(length1,1) -        href(1,1);
dx2 = idealPower2(length1,1) - idealPower2(1,1);
m2= dy2/dx2;
root2 = -b2/m2