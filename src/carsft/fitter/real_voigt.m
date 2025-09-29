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