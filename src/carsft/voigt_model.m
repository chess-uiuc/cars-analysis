function [chi,area,alphad] = voigt_model(chi1,wraman,waven,gamma,T,cmass)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function calculated the Raman-resonant susceptibility using a Voigt
%lineshape to account for collisional dephasing without velocity changing
%collisions and Doppler broadening
%
% chi = susceptibility vector for a single gas (unscaled)
%
% chi1 = prepopulated transition amplitude = (cross section) X (pop diff)
% wraman = Raman transition frequencies in cm-1
% waven = wavenumber grid for susceptibility calculation (cm-1)
% gamma = collision widths in cm-1
% T = temperature (K)
% cmass = mass of gas species in amu (Daltons)


c = 3e10; %vacuum speed of light in cm-1
kb = 1.38e-23; %Boltzmann constant in J/K
mo = 1.660539e-27; %convert 1 amu to kg (Daltons to kg)

mass = cmass*mo; %mass of a single molecule in kg
alphad0 = (1/c)*sqrt(2*kb*T/mass)*100; 
chi = zeros(length(chi1),length(waven));

for k = 1:length(chi1)  %loops through Raman transitions (better with bsxfun ?)
    alphad = alphad0*wraman(k); %Doppler 'sigma' width from CARSFT (in cm-1)
    xv = (wraman(k)-waven)/alphad;
    yv = gamma(k)/(2*alphad);
    z = xv+i*yv;
    F = faddeeva(z,16);

    ac = sqrt(pi)/(2*alphad);
    chi(k,:) = chi1(k)*ac*(imag(F) + i*real(F));
    %area = trapz((waven-wraman(k)),imag(chi(k,:)))/(sqrt(pi)*alphad);
end

chi = sum(chi,1);

return