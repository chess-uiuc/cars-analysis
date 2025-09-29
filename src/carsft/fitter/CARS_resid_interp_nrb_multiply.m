function f = CARS_resid_interp_nrb_multiply(x);

%%%%This script computes the residual for a 3-parameter (T, CO, XNR) fit to
%%%%CARS spectra using calls to full frequency-domain calculation
%%%%This version multiplies the theory by the measured argon NRB spectrum
global data omega_grid NRB
global wlib CARS Tlib Plib Alib

T=x(1);
CO=x(2);
alpha=x(3);
hshift=x(4);
wexp=x(5);
vshift=x(6);

dtp = 8000;dtau = 0;dwp = 1;P = 1;

theory = interpolate3d_lagrange(CARS,wlib,Tlib,Plib,Alib,T,CO,alpha,wexp*(omega_grid-hshift));
%theory = freq_domain_cars_multi_v2(wexp*(omega_grid-hshift),T,P,[1 CO],dtp,dtau,alpha,dwp);
%theory=interpolate2d(CARS,wlib,Tlib,Plib,T,P,wexp*(omega_grid-hshift));
%theory=interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,P,wexp*(omega_grid-hshift));

theory=theory.*NRB;
theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));
%theory=theory/max(theory);

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end

indx=find(data ~= 0);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));  %CHANGE
%THIS BACK FOR REAL DATA!!
%f=((data(indx)-vshift)-theory(indx))./sqrt(abs(theory(indx)+1e-5));
