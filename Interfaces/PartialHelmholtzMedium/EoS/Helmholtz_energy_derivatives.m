%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives of the Helmholtz energy
% Matthis Thorade, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
syms pi1 pi2 pi3 bi1 bi2 bi3 bi4 gi1 gi2 gi3 gi4 gi5 gi6 gi7 gi8 gi9;

% Terms for Helmholtz energy
logTerms = + li1*log(tau^li2);
idealPolyTerms = + pi1*tau^pi2;
EinsteinTerms = + ei1*log(1 - exp(ei2*tau));
coshTerms = - ci1*log((cosh(ci2*tau)));
sinhTerms = + si1*log((sinh(si2*tau)));
% fi = +log(delta) +logTerms +idealPolyTerms +EinsteinTerms +coshTerms +sinhTerms;
fi = +sinhTerms;

residualPolyTerms = pi1*tau^pi2*delta^pi3;
BWRTerms = bi1*tau^bi2*delta^bi3*exp(-delta^bi4);
GaussTerms = gi1*tau^gi2*delta^gi3*exp(gi6*(delta-gi9)^2 + gi7*(tau-gi8)^2);
nonAnalyticalTerms = 0;
% fr = +residualPolyTerms +BWRTerms +GaussTerms +nonAnalyticalTerms;
fr = +GaussTerms; 


% First derivatives
fit = diff(fi,tau);
fit = simplify(fit)
fid = diff(fi,delta);
fid = simplify(fid);

frt = diff(fr,tau);
frt = simplify(frt);
frd = diff(fr,delta);
frd = simplify(frd);


% Second derivatives
fitt = diff(fit,tau);
fitt = simplify(fitt)

frtt = diff(frt,tau);
frtt = simplify(frtt);
frtd = diff(frt,delta);
frtd = simplify(frtd);
frdd = diff(frd,delta);
frdd = simplify(frdd);


% Third derivatives
fittt = diff(fitt,tau);
fittt = simplify(fittt)

frttt = diff(frtt,tau);
frttt = simplify(frttt)
frttd = diff(frtd,tau);
frttd = simplify(frttd);
frtdd = diff(frtd,delta);
frtdd = simplify(frtdd);
frddd = diff(frdd,delta);
frddd = simplify(frddd);
