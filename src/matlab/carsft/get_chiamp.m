function [chiamp] = get_chiamp_v2(gasname,molpar,T,Tv,X,QN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute amplitudes of CHI1111 (this code is currently valid for
% parallel-polarized laser beams only
% INPUTS:
% 1. gasname = names of gases with nonzero mole fraction from cars.mol
% 2. molpar = structure array of molecular parameters from cars.mol
% 3. T = rotational temperature (K)
% 4. Tv = vibrational temperature of species 'gas' (K)
% 5. X = mole fraction vector for species with user input X > 0
% 6. QN = structure array from get_trans.m: contains Qo = [vmax Jmax] vi Ji
%           branch names
%
% OUTPUTS:
% chiamp = vector of transition amplitudes CHI1111 (length = Ntrans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chiamp = []; %initiallize the output variable

for k = 1:length(gasname)
  kgas = find(molpar.gasname == gasname(k));

  gI = molpar.gI(kgas,:);we = molpar.WE(kgas);Be = molpar.BE(kgas);
  a2 = molpar.AC1(kgas); g2 = molpar.AG(kgas)*a2;
  a2 = a2*a2;g2 = g2*g2; %mean and anisotropy of the polarizability derivate (for Q,S,O transitions)
  gam2 = molpar.GAM(kgas);gam2 = gam2.*gam2; %anisotropy of the mean polarizability (for ROT transitions)

  qo = QN{k}.Qo;vi = QN{k}.vi;Ji = QN{k}.Ji; %initial (v,J) quantum numbers for Raman transitions
  branch = QN{k}.branch;

  CHIAMP = zeros(1,length(QN{k}.Ji));

  bJJQ=Ji.*(Ji+1)./((2*Ji-1).*(2*Ji+3));  %Placzek-Teller coeffcients, Q branch
  bJJS=1.5*(Ji+1).*(Ji+2)./((2*Ji+1).*(2*Ji+3)); %S branch and ROT
  bJJO=1.5*Ji.*(Ji-1)./((2*Ji+1).*(2*Ji-1)); %O branch 

  %1. compute the required Boltzmann fractions up to J = Jmax +2, v = vmax +1
  NJo = J_fraction([0:1:qo(2)+2],Be,gI,T);
  Nvo = v_fraction(qo(1)+1,we,Tv(kgas));

  %2. compute the delta-N terms for each branch
  indxQ = find(strcmp(branch,'Q') == 1);
  if (isempty(indxQ) ~= 1)
      nj1 = NJo(Ji(indxQ)+1);nv1 = Nvo(vi(indxQ)+1);nv2 = Nvo(vi(indxQ)+2);
      DNQ = nj1.*(nv1 - nv2);
      CHIAMP(indxQ) = DNQ.*(vi(indxQ)+1).*(a2+g2*bJJQ(indxQ)*4/45);
  end

indxS = find(strcmp(branch,'S') == 1);
if (isempty(indxS) ~= 1)
    nj1 = NJo(Ji(indxS)+1);nj2 = NJo(Ji(indxS)+3);  
    nv1 = Nvo(vi(indxS)+1);nv2 = Nvo(vi(indxS)+2);
    DNS = nv1.*nj1 - nv2.*nj2;
    CHIAMP(indxS) = DNS.*(vi(indxS)+1).*(g2*bJJS(indxS)*4/45);
end

indxO = find(strcmp(branch,'O') == 1);
if (isempty(indxO) ~= 1)
    nj1 = NJo(Ji(indxO)+1);nj2 = NJo(Ji(indxO)-1);
    nv1 = Nvo(vi(indxO)+1);nv2 = Nvo(vi(indxO)+2);
    DNO = nv1.*nj1 - nv2.*nj2;
    CHIAMP(indxO) = DNO.*(vi(indxO)+1).*(g2*bJJO(indxO)*4/45);
end

indxROT = find(strcmp(branch,'ROT') ==1);
if (isempty(indxROT) ~= 1)
    nj1 = NJo(Ji(indxROT)+1);nj2 = NJo(Ji(indxROT)+3);
    nv1 = Nvo(vi(indxROT)+1);
    jfac = (2*Ji(indxROT)+1)./(2*(Ji(indxROT)+2)+1);
    DNROT = nv1.*(nj1 - jfac.*nj2);
    CHIAMP(indxROT) = DNROT.*(gam2*bJJS(indxROT)*4/45);
end

chiamp = [chiamp; X(k)*CHIAMP'];

end
