%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,theory] = CARS_resid_2Dlib_interp(x);
%%%%This script computes the residual for a 2-parameter (e.g. T, XNR) fit to
%%%%CARS spectra using calls to the LaGrange interpolated 2D library
%%%%This version multiplies the theory by the measured argon NRB spectrum
%%%%Simply make NRB = ones(1,length(omega_grid)) to fit to argon-normalized
%%%%data

global data omega_grid NRB
global wlib CARS Tlib Plib Alib
global weights

T=x(1);
v2=x(2);
hshift=x(3);
wexp=x(4);
vshift=x(5);

theory = interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,v2,wexp*(omega_grid-hshift));
theory=theory.*NRB;
theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end
indx=find(data ~= -1e5);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-10));

f = weights.*f;

end