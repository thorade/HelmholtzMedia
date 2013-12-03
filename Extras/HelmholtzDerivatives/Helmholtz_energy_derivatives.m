%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives of the Helmholtz energy
% Matthis Thorade, 2012%
% warning: this might keep MatLab busy for a while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% Helmholtz energy
syms fi fr;
% Helmholtz energy first derivatives
syms fit frt fid frd; 
% Helmholtz energy second derivatives
syms fitt frtt fidd frdd fidt frdt;
% Helmholtz energy third derivatives
syms fittt frttt fiddd frddd;

% independent variables
syms delta tau;

% Terms
syms logTerms idealPolyTerms EinsteinTerms coshTerms sinhTerms;
syms residualPolyTerms BWRTerms GaussTerms nonAnalyticalTerms;

% Parameters in the bank of terms
syms li1 li2 ei1 ei2 ci1 ci2 si1 si2;
syms pi1 pi2 pi3;
syms bi1 bi2 bi3 bi4;
syms gi1 gi2 gi3 gi4 gi5 gi6 gi7 gi8 gi9;
syms ni1 ni2 ni3 ni4 ni5 ni6 ni7 ni8 ni9 ni10 ni11 ni12 Distance Phi;

% Terms for Helmholtz energy
logTerms = + li1*log(tau^li2);
idealPolyTerms = + pi1*tau^pi2;
EinsteinTerms = + ei1*log(1 - exp(ei2*tau));
coshTerms = - ci1*log((cosh(ci2*tau)));
sinhTerms = + si1*log((sinh(si2*tau)));
% fi = +log(delta) +logTerms +idealPolyTerms +EinsteinTerms +coshTerms +sinhTerms;
fi = +EinsteinTerms ;
fi = simplify(fi)

residualPolyTerms = pi1*tau^pi2*delta^pi3;
BWRTerms = bi1*tau^bi2*delta^bi3*exp(-delta^bi4);
GaussTerms = gi1*tau^gi2*delta^gi3*exp(gi6*(delta-gi9)^2 + gi7*(tau-gi8)^2);
Distance = ((1-tau)+ni8*((delta-1)^2)^(1/(2*ni7)))^2 + ni11*((delta-1)^2)^ni12
Phi = exp(-ni9*(delta-1)^2-ni10*(tau-1)^2)
nonAnalyticalTerms = ni1*delta*Distance^ni6*Phi;
% fr = +residualPolyTerms +BWRTerms +GaussTerms +nonAnalyticalTerms;
fr = +GaussTerms  
fr = simplify(fr)

%% First derivatives
fit = diff(fi,tau);
fit = simplify(fit);
fid = diff(fi,delta);
fid = simplify(fid);

frt = diff(fr,tau);
frt = simplify(frt)
frd = diff(fr,delta);
frd = simplify(frd)

%% Second derivatives
fitt = diff(fit,tau);
fitt = simplify(fitt);

frtt = diff(frt,tau);
frtt = simplify(frtt)
frtd = diff(frt,delta);
frtd = simplify(frtd)
frdd = diff(frd,delta);
frdd = simplify(frdd)

%% Third derivatives
fittt = diff(fitt,tau);
fittt = simplify(fittt);

frttt = diff(frtt,tau);
frttt = simplify(frttt)
frttd = diff(frtd,tau);
frttd = simplify(frttd)
frtdd = diff(frtd,delta);
frtdd = simplify(frtdd)
frddd = diff(frdd,delta);
frddd = simplify(frddd)

%% ideal gas limit
%limit(frd,delta,0,'right')
