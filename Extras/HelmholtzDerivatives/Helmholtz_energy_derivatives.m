%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives of the Helmholtz energy EoS
% Matthis Thorade, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

% symbolic variables for Helmholtz energy
syms fi fr delta tau;
% symbolic variables for bank of terms
syms li1 li2 ei1 ei2 ci1 ci2 si1 si2;
syms pi1 pi2 pi3;
syms bi1 bi2 bi3 bi4;
syms gi1 gi2 gi3 gi4 gi5 gi6 gi7 gi8 gi9;
syms ni1 ni2 ni3 ni4 ni5 ni6 ni7 ni8 ni9 ni10 ni11 ni12 Distance Phi;
% terms for Helmholtz energy
logTerms = + li1*log(tau^li2);
idealPolyTerms = + pi1*tau^pi2;
EinsteinTerms = + ei1*log(1 - exp(ei2*tau));
coshTerms = - ci1*log((cosh(ci2*tau)));
sinhTerms = + si1*log((sinh(si2*tau)));
residualPolyTerms = pi1*tau^pi2*delta^pi3;
BWRTerms = bi1*tau^bi2*delta^bi3*exp(-delta^bi4);
GaussTerms = gi1*tau^gi2*delta^gi3*exp(gi6*(delta-gi9)^2 + gi7*(tau-gi8)^2);
Distance = ((1-tau)+ni8*((delta-1)^2)^(1/(2*ni7)))^2 + ni11*((delta-1)^2)^ni12
Phi = exp(-ni9*(delta-1)^2-ni10*(tau-1)^2)
nonAnalyticalTerms = ni1*delta*Distance^ni6*Phi;


%% Helmholtz energy
% fi = +log(delta) +logTerms +idealPolyTerms +EinsteinTerms +coshTerms +sinhTerms;
fi = +EinsteinTerms ;
fi = simplify(fi)
% fr = +residualPolyTerms +BWRTerms +GaussTerms +nonAnalyticalTerms;
fr = +BWRTerms  
fr = simplify(fr)

%% Helmholtz energy first derivatives
fit = diff(fi,tau);
fit = simplify(fit);
fid = diff(fi,delta);
fid = simplify(fid);
frt = diff(fr,tau);
frt = simplify(frt)
frd = diff(fr,delta);
frd = simplify(frd)

%% Helmholtz energy second derivatives
fitt = diff(fit,tau);
fitt = simplify(fitt);
frtt = diff(frt,tau);
frtt = simplify(frtt)
frtd = diff(frt,delta);
frtd = simplify(frtd)
frdd = diff(frd,delta);
frdd = simplify(frdd)

%% Helmholtz energy third derivatives
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
