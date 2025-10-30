function [wraman,gamma,qn] = get_trans_v3(gasname,molpar,T,Tv,wexp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
% 1. gasname = names of gases in the calculation, obtained from cars.mol
% 2. molpar = cell array containing parameters from cars.mol
% 3. T = rotational temperature (K)
% 4. Tv = vibrational temperature for gas (K)
% 5. wexp = experimental wavenumber grid (cm-1)
%
% Outputs
% 1. wraman = frequencies of thermally populated lines within
%   +/-N*gamma(v,J) of the limits of the experimental wavenumber grid
% 2. gamma = Raman collisional linewidths within +/-N*gamma(v,J)....
% 3. qn = quantum numbers and Raman branches contained in structure array
%    qn.Qo = [vmax Jmax] max quantum numbers retained in the calc. for gas
%    qn.vi = initial v quantum number for returned transitions Q,S,O,ROT
%    qn.Ji = inital J quantum number for transition Q,S,O,ROT
%    qn.branch = string varible names of Raman branch for returned transitions
%              Q.S,O,ROT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 50; %include transitions +/-N*gamma(v,J) from experimental wave# limits

wraman = []; gamma = []; qn = {}; %initialize output variables

for k = 1:length(gasname)
  kgas = find(molpar.gasname == gasname(k));
  Be=molpar.BE(kgas);we=molpar.WE(kgas);gI=molpar.gI(kgas,:);

  %A. calculate Boltzmann dist. for 99.5% of the thermal population
  vo = [0:1:15];Jo = [0:1:200];
  NJo = J_fraction(Jo,Be,gI,T);Nvo = v_fraction(max(vo),we,Tv(kgas));
  [Nv,vmax] = Boltzmann_limit(Nvo);[NJ,Jmax] = Boltzmann_limit(NJo);
  vo = [0:1:vmax];Jo = [0:1:Jmax];QN.Qo = [vmax Jmax];

  formatSpec = 'maximum vibrational quantum number for %s = %2u\n';
  fprintf(formatSpec,gasname(k),vmax);
  formatSpec = 'maximum rotational quantum number for %s = %2u\n';
  fprintf(formatSpec,gasname(k),Jmax);
  fprintf('\n')

  %B. calculate the raman shift for transitions at vi,Ji accounting for
  % 99.5% of the thermal population
  [wQ,wS,wO,wROT] = raman_shift(gasname(k),molpar,vmax,Jmax);
  WRAMAN = [wQ' wS' wO' wROT'];

  %C. %%%compute the Raman linewidths in FWHM cm-1
  [gammaQ gammaS gammaO gammaROT] = meg_linewidths(Jo,T,gasname(k));
  gammaQ = repmat(gammaQ,1,vmax+1);gammaS = repmat(gammaS,1,vmax+1);
  gammaO = repmat(gammaO,1,vmax+1);gammaROT = repmat(gammaROT,1,vmax+1);
  GAMMA = [gammaQ gammaS gammaO gammaROT];

  %D. assign branch names each transition
  branch = [repmat({'Q'},Jmax+1,vmax+1)...
                repmat({'S'},Jmax+1,vmax+1)...
                repmat({'O'},Jmax+1,vmax+1)...
                repmat({'ROT'},Jmax+1,vmax+1)];

  %E. assign initial quantum numbers to each transition
  vi = [0:1:vmax];vi = repmat(vi,Jmax+1,1);vi = repmat(vi,1,4);
  Ji = transpose([0:1:Jmax]);Ji = repmat(Ji,1,vmax+1);Ji = repmat(Ji,1,4);
  indx = find(strcmp(branch,'O') == 1);Ji(indx) = Ji(indx) + 2;  %shift Ji for O-branch by +2


  %E. %%%determine which transitions to keep based on 
  % +-N*gamma from the wavenumber limits of the calculation
  indx = find( (WRAMAN - N*GAMMA) <= max(wexp) & (WRAMAN + N*GAMMA) >= min(wexp) );
  WRAMAN = WRAMAN(indx);GAMMA = GAMMA(indx);
  vi = vi(indx);Ji = Ji(indx);gas = repmat(gasname(k),length(indx),1); branch = branch(indx);


  wraman = [wraman; WRAMAN]; gamma = [gamma; GAMMA];
  qn{k}.vi = vi; qn{k}.Ji = Ji; qn{k}.gas = gas; qn{k}.Qo = [vmax Jmax]; qn{k}.branch = branch;

end

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

if (gas == 'H2')
    gamma = meg_co(J,T);
    gamma = zeros(size(gamma))+0.01;
end

%%%compute S- and O-branch widths from RW approximation
gammaS = 0.5*(gamma + circshift(gamma,-2,1));
gammaS = gammaS(1:length(J)-2);

gammaO = gammaS; %O-branch values (note: Ji vector for O starts at J = 2)

gamma = gamma(1:length(J)-2); %Q-branch values

gammaROT = gammaS; %pure-rotational S-branch values

end