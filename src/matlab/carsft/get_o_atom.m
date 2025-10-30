function [chiamp,gamma,wraman,N,DN] = get_o_atom(T,molpar)

%function returns the frequency-dependent susceptibility for O atom
%considers only 3P2-3P0 and 3P2-3P1 transitions. Dasch and Bechtel did not
%observe 3P1-3P0 in their experiments and estimate its cross section to be
%small

%calculate/estimate the polarizability using results from Dasch and Bechtel
%(1981). Numbers in their Table 1 are relative to pure-rotational S(7)
%transition in O2.

o2col = find(strcmp(molpar.gasname,"O2") == 1); %locate molecular O2 data from cars.mol
gam2 = molpar.GAM(o2col);gam2 = gam2*gam2; %polarizability anisotropy for O2

Jref = 7; %reference S(7) line in Dasch and Bectel
gam_ref = gam2*(2*Jref+1)/(2*(Jref+2)+1)*(4/45); %linestrength of O2 S(7) in polarizability units
sigma_ref = 1.4e-29; %cross section for S(7) ref. transition in cm2/sr from Dasch and Bechtel

sigma(1) = 6e-31;sigma(2) = 4e-31; %relative sigma_zz for 3P2-3P1 and 3P2-3P0 transitions
sigma(1) = gam_ref*sigma(1)/sigma_ref; %linestrengths for O atom transitions in polarizability units
sigma(2) = gam_ref*sigma(2)/sigma_ref; %these are the units used in the main code

N = o_atom_fraction(T); %Boltzmann population of O atom electronic states
DN(1) = N(1)-(5/3)*N(2);
DN(2) = N(1)-(5)*N(3); %multiply by (2J+1)/(2J'+1) factors (check this for O atom?)

chiamp = DN.*sigma; %amplitude of chi

gamma = o_atom_relax;

wraman = [158.265 226.977];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = o_atom_fraction(T);

M = 5; %number of states in the 1s2-2s2-2p4 ground configuration of O atom
        %1 = 3P2, degeneracy = 2(2)+1 = 5, energy = 0
        %2 = 3P1, degeneracy = 2(1)+1 = 3, energy = 158.265 cm-1
        %3 = 3P0, degeneracy = 2(0)+1 = 1, energy = 226.977 cm-1
        %4 = 1D2, degeneracy = 2(2)+1 = 5, energy = 15,867.862 cm-1
        %5 = 1S0, degeneracy = 2(0)+1 = 1, energy = 33,972.582 cm-1
kb = 0.69503476; %Boltzmann constant in cm-1/K        

J = [2 1 0 2 0]; %total (spin + orbital) angular momenta
energy = [0 158.265 226.977 15867.862 33972.582]; %level energies
        
q = (2*J+1).*exp(-energy/(kb*T));

N = q/sum(q);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = o_atom_relax

%returns linewidths based on O-atom relaxation data (or other information)

%relaxation times "tau_CARS" in ps, measured in Ar/O2 plasma at Sandia/NM
%at P = 0.82 atm
tau(1) = 139;
tau(2) = 127;

%convert to linewidths in cm-1/atm (FWHM)
co = 2.9999e10; %speed of light in cm/s
gamma = 1./(2*pi*co*tau*1e-12); %linewidth cm-1 FWHM;
gamma = gamma/0.82; %linewidth cm-1/atm FWHM;

end

