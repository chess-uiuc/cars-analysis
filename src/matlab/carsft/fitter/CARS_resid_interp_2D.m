function [f,theory] = CARS_resid_interp_2D(x);
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


%1. compute the slit function from the x(6) and x(7) inputs
if (length(x) == 6)  %case 1: Gaussian slit
    dwG = x(6); %Gaussian slit width FWHM (cm-1)
    maxws = 10*dwL;dws = mean(diff(wlib));ws = [-maxws:dws:maxws]; %wave# grid for slit convolution
    slit = exp(-4*log(2)*(ws/dwL).^2);
end

if (length(x) == 7); %case 2: voigt slit
    dwG = x(6); %Gaussian slit contribution FWHM (cm-1)
    dwL = x(7); %Lorentz slit contribution FWHM (cm-1)
    maxws = 10*max(dwG,dwL);dws = mean(diff(wlib));ws = [-maxws:dws:maxws]; %wave# grid for slit convolution
    slit = real_voigt(ws,dwG,dwL); 
end

%2. limit wlib to the bounds of omega_grid. This let's us apply wexp and
%hshift
wtheory = [min(omega_grid):dws:max(omega_grid)];

%3. interpolate the theoretical susceptibility from the library
%theory = interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,v2,wexp*(omega_grid-hshift));
theory = interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,v2,wtheory);
theory = conv_spk(slit,theory.*conj(theory));  %NOTE: slit fxn is computed once on unshifted/unstretched grid. this is ok for small wexp 
theory = interp1(wtheory,theory,wexp*(omega_grid-hshift));
theory=theory.*NRB;
theory=theory/abs(trapz(omega_grid,theory));
theory(find(isnan(theory)==1))=0;

%theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end
indx=find(data ~= 0);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));
f = weights.*f;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = real_voigt(x,sigma,gamma)

%function to compute the real-valued voigt function on 1D grid, x
%
% x = wave# grid
% sigma = gaussian width parameter (std. dev)
% gamma = lorentz width parameter (FWHM)

xv = x/sigma;
yv = gamma/(2*sigma);
z = xv+i*yv;
F = faddeeva(z,16);

f = real(F);

end