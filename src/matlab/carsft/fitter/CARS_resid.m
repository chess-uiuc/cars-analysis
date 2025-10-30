function f = CARS_resid(x);

%%%%This script computes the residual for a 3-parameter (T, CO, XNR) fit to
%%%%CARS spectra using calls to full frequency-domain calculation
global data omega_grid NRB
global wlib CARS Tlib Plib

T=x(1);
CO=x(2);
alpha=x(3);
hshift=x(4);
wexp=x(5);
vshift=x(6);

dtp = 8000;dtau = 0;dwp = 1;P = 1;

theory = interpolate3d_lagrange(CARS,wlib,Tlib,
%theory = freq_domain_cars_multi_v2(wexp*(omega_grid-hshift),T,P,[1 CO],dtp,dtau,alpha,dwp);
%theory=interpolate2d(CARS,wlib,Tlib,Plib,T,P,wexp*(omega_grid-hshift));
%theory=interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,P,wexp*(omega_grid-hshift));

theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));
%theory=theory.*NRB;
%theory=theory/max(theory);

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end

indx=find(data ~= 0);
%f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));  %CHANGE
%THIS BACK FOR REAL DATA!!

f=((data(indx)-vshift)-theory(indx))./sqrt(abs(theory(indx)+1e-5));
