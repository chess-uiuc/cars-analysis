function f = CARS_resid_3Dlib_interp(x);

%%%%This script computes the residual for a 2D parameterized library of
%%%%CARS spectra
global data omega_grid NRB
global wlib CARS Tlib Plib,Xlib

T=x(1);
P=x(2);
X = x(3);
hshift=x(4);
wexp=x(5);
vshift=x(6);

theory=interpolate3d_lagrange(CARS,wlib,Tlib,Plib,T,P,wexp*(omega_grid-hshift));
theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));
%theory=theory.*NRB;
%theory=theory/max(theory);

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end

indx=find(data ~= 0);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));
