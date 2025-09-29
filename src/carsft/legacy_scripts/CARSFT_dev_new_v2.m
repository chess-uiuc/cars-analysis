function [S,CHI,w,E3,w3] = CARSFT_dev_new_v2(wexp,T,P,X,dtp,dtau3,alpha,dwp) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency domain CARS calculation for femtosecond/impulsive Raman preparation, 
% with isolated Raman lines
%
% wexp = "experimental" wavenumber grid (cm-1); output is interpolated onto
%        wexp
% T = temperature (K) (if scalar then equilibrium, if vector then...
%        T(1) = Trot
%        T(2) = Tv for gas 1; T(3) = Tv for gas 3 etc...
% P = pressure (atm)
% X = vector of mole fractions [N2 gas2 gas3 gas 4] 
% dtp = linewidth of Gaussian probe laser (cm-1)
% dtau3 = time delay of probe laser (ps)
% alpha = scaling factor for nonresonant suceptibility Xnr = alpha*(Xr(v =
% 1)
% dwp = slit width (assumed Gassian, in cm-1)
%
% S = CARS spectrum, interpolated onto input wavenumber grid, wexp
% CHI = frequency-domain CARS susceptibility (not convolved with probe)
% E3 = probe pulse field envelope*exp(i*w*dtau3) (computed from transform
% limit of dtp)
% w3 = wavenumber grid (cm-1) for probe pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%computational wavenumber grid
wmin = 2000;wmax = 2400;dw = 0.005;  %hard wire the wavenumber grid for now
                                    %wavenumber step = 0.01 cm-1 << Raman
                                    %linewidth
w = [wmin:dw:wmax];    %computational frequency grid in cm-1
CHI = zeros(1,length(w)); %dedicate memory for susceptibility vector

%%%some constants
kb = 1.38e-23; %Boltzmann constant (J/K)
c = 2.99999e10; %speed of light, vacuum (cm/s)
mo = 1.660539e-27; %convert 1 amu to kg (Daltons to kg)
voigt_threshold = 1; %thereshold ratio of Doppler to collision width that throws flag for voigt lineshape model

%%%compute the "effective" probe field
dw3 = 14.7/dtp; %probe linewidth in cm-1
w3 = [-20:dw:20];
E3 = exp(-2*log(2)*(w3/dw3).^2);E3 = E3.*exp(i*w3*dtau3*2*pi*3e10/1e12);
%Note the sign in the argument of exp above is flipped. Marrocco's
%expression for Yuratich Eq. show -i*w3... This could be to compensate for
%a factor of i or a sign flip in CHI below. In any case this code agrees
%with the direct time-domain calculation at T = 2000 K, P = 1 atm, for 50,
%100, and 150 ps FWHM probe pulses at delays up to 6X the probe pulse width
%using a wavenumber grid spacing of dw = 0.01 cm-1. There are some sligt
%differences comparing to time-domain calculation for a 6-ps probe.

%%read gas parameters from CARS.MOL
[molpar] = read_cars_dot_mol('CARS.mol');

gasindx = find(X ~= 0); %find the gas species to include in the calculation
gasname = molpar.gasname(gasindx);

%assign rot and vib temperatures based on structure of T variable
[T,Tv] = assign_temps(T,length(molpar.gasname));
gamma = [];wraman = []; qn = {}; chiamp = [];


for kk = 1:length(gasindx)  %loop through the gas species to be considered

kgas = gasindx(kk);

%get the raman frequencies, linewidths, and branches of lines to be
%included in the calculation: (1) thermally populated; (2) within
%N*gamma(J) of wexp limits
[wraman,gamma,branch,qn] = get_trans_v2(gasname,molpar,T,Tv,wexp);   

%compute the Raman resonant amplitude of each line
[chiamp] = get_chiamp_v2(gasname,molpar,branch,T,Tv,qn);

end


%%%representative Doppler width in cm-1 based on average transition
%%%frequency
alphad = (1/c)*sqrt(2*kb*T/(cmass*mo))*100; %this is not the FWHM, its the 'ALPHAD' variable in CARSFT
alphad = alphad*mean(wraman);  

%%%comapare Doppler width estimate to collision widths and select lineshape model
ratiodplr = alphad./gamma;ratiodplr = max(ratiodplr);

if (ratiodplr < voigt_threshold); linemodel = 'iso';
else; linemodel = 'voigt';end

formatSpec = 'lineshape model for %s = %s\n';
fprintf(formatSpec,gas(kgas),linemodel);
fprintf(' \n');

if (isequal(linemodel,'iso') == 1);  %default model is isolated lines in main code
    for k = 1:length(DN)
        waven = wraman(k)-w;
        chi = chi + chi1(k)./(waven - i*gamma(k)/2);
    end
end

if (isequal(linemodel,'voigt') == 1);
    chi = voigt_model(chi1,wraman,w,gamma,T,cmass);
end

chi = chi+alpha*dsdw(1);  %add the nonresonant component (you need to change this. the way it is...
                            %currently written, the nonresonant background
                            %is referenced independently as a fraction of
                            %all species. We'd lilke be better of
                            %referencing to N2 (or a 'ref' gas)

CHI = CHI+X(kgas)*chi;    %update the TOTAL gas mixture susceptibility, CHI 

end %end loop through gas species

CHI = (P/T)*CHI;


%%%convolve CARS susceptibility with probe pulse
wd = mean(diff(w));w3a=[min(w3):wd:max(w3)]; %wavenumber grid for probe and slit convolutions
E3 = interp1(w3,E3,w3a);  %interpolate probe onto wd grid
S = conv(CHI,E3,'same');S = S.*conj(S);
%S = conv_spk(CHI,E3);S = S.*conj(S);

%%%convolve CARS spectrum with slit function
ifunc=exp(-4*log(2)*(w3a/dwp).^2);  %instrument function (assumed Gaussian)
S = conv(S,ifunc,'same');
%S = conv_spk(S,ifunc);

%%%interpolate the output onto the "experimental" wavenumber grid, wexp
S = interp1(w,S,wexp);S = S/max(S);

end

%%%%%%%%%%%%%%%%%%%addtional subroutines%%%%%%%%%%%%%%%%%%%
function [gamma,gammaS,gammaO,gammaROT] = meg_linewidths(J,T,gas);

J = [J max(J)+1 max(J+2)]; %expand J vector to compute S- and O-branch from RWA

if (gas == 'N2')
    gamma=meg_n2_o2(J,T);gamma=gamma(:,1);
end

if (gas == 'O2')
    gamma=meg_n2_o2(J,T);gamma=gamma(:,2);
end

if (gas == 'CO')
    gamma = meg_co(J,T);
end

%%%compute S- and O-branch widths from RW approximation
gammaS = 0.5*(gamma + circshift(gamma,-2,1));
gammaS = gammaS(1:length(J)-2);

gammaO = 0.5*(gamma + circshift(gamma,2,1));
gammaO = [0;0;gammaO(3:length(J)-2)];

gamma = gamma(1:length(J)-2); %Q-branch values

gammaROT = gammaS; %pure-rotational S-branch values

end

function [dsdw] = raman_cross_section(gas,v,J,wp,eQ,M,a2,g2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computes the Raman cross section FOR Q BRANCH transition using Eqs. 3 and 4 in
%Lucht, Palmer, and Maris, Opt. Lett. 12, 386 (1987)
%
%dsdw = (v+1)*(a^2 + (4/45)*bJJ*gamma^2)
%
%INPUTS:
%  1. AG = d(gamma)/dQ / da/dQ = ratio of polarizability derivate anisotropy to
%       mean polarizabilty derivative
%  2. AC1 = mean polarizability derivative (da/dQ)/sqrt(2*pi*c*wRAMAN*reduced_mass)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = zeros(1,lenght(eQ));
c = 3e10; %speed of light
nv = length(v);

%%%reduced mass calculation
if (gas == 'N2')
    mu = 1/14.0067 + 1/14.0067; mu = 1/mu;
end

if (gas == 'CO')
    mu = 1/12.0107 + 1/15.999; mu = 1/mu;
end

if (gas == 'O2')
    mu = 1/15.999 + 1/15.999; mu = 1/mu;
end

if gas == 'H2'
    mu = 1/1.00784 + 1/1.00784; mu = 1/mu;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = ac1*sqrt(mu*eQ); %this includes the 1/wRaman dependence
a = ac1*sqrt(mu); %multiply by reduced mass correction does CARSFT do this?
g = ag*a; 
a = a.*a;g = g.*g;

a=1.91e-29;  %mean polarizability derivative for N2 (cm4)
g=7.10e-30;


bJJ = J.*(J+1)./((2*J-1).*(2*J+3));bJJ=transpose(bJJ); 
bJJ = repmat(bJJ,length(v)-1,1);

sigma = a + (4/45)*bJJ.*g;
sigma = ((wL+eQ).^4).*sigma;
sigma = reshape(sigma,length(J),length(v)-1);
sigma = bsxfun(@times,v(1:nv-1)+1,sigma);
sigma = reshape(sigma,nJ*(nv-1),1);

end

function [N,indx] = Boltzmann_limit(No)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function computes the max v and max J for the calculation 
%based on threshold total population
%
% No = Boltzmann distribution for large set of quantum numbers
%       greater than what would be needed for temperatures of interest
%
% N = truncated Boltzmann distribution to contain threshold level of 
%       population
% indx = maximum quantum number in truncated distribution, such that 
%         quantum number = [0:1:indx] is size indx+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = 0.995; %threshold for fraction of total population retained

indx = find(cumsum(No) <= threshold);

if (max(indx) < length(No))
    indx = max(indx);
else                    
    indx = max(indx);
end
    
if (length(indx) < 1)
    indx = 2;
end

N = No(1:indx+1);

end


function chi = voigt_model(chi1,wraman,waven,gamma,T,cmass)

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

    chi(k,:) = chi1(k)*(imag(F) + i*real(F));
    chi(k,:) = chi(k,:)/(sqrt(pi)*alphad);
end

chi = sum(chi,1);

end

function [Tr,Tv] = assign_temps(T,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function assigns rotational and vibrational temperatures from the
% input temperature array
% 
% T = input array from user (K)
% N = length(gas) = number of species in the calculation
% Tr = single rotational temperature (K)
% Tv = vibrational temperatures for all gases in cars.mol (K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tr = T(1); %single rotational temperature
Tv = zeros(N,1); %vibrational temperatures for each gas in CARS.mol
if length(T) == 1; Tv=Tv+Tr;else
    Tv(1:length(T)-1) = T(2:length(T));
    indx = find(Tv==0);Tv(indx)=Tr;
end

if (length(T)>1 & length(T)~=N+1)
    formatSpec = '\n warning: size of T array does not match no. of gases in CARS.mol \n';
    fprintf(formatSpec);
    formatSpec = ' assigning remaining gases to Trot \n';
    fprintf(formatSpec);fprintf('\n');
end

end



