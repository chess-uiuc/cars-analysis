function [omegaQ,omegaS,omegaO,omegaROT] = raman_shift(gas,molpar,maxV,maxJ)

v=[0:1:maxV+1];J=[0:1:maxJ+2];

k = find(molpar.gasname == gas);
we=molpar.WE(k);wexe=molpar.WX(k);weye=molpar.WY(k);weze=molpar.WZ(k);
Be=molpar.BE(k);alphae = molpar.ALPHAE(k);gammae = molpar.GAME(k);
De=molpar.DE(k);betae=molpar.BETAE(k);deltae=molpar.DELTE(k);
Ho=molpar.H0(k);He=molpar.HE(k);

wo = we-wexe+0.75*weye;woxo = wexe-1.5*weye;woyo = weye;
gov = wo*v-woxo*v.^2+woyo*v.^3;

%%%%calculates v- and J-dependent rotational term energies, FJ%%%%
VPH=(v+0.5);RJ=J.*(J+1);
C1 = Be+VPH.*(-alphae+VPH*gammae);
C2 = -(De+VPH.*(betae+VPH*deltae));
C3 = Ho+VPH*He;

FJ=zeros(length(v),length(J));GvJ=zeros(length(v),length(J));
for k = 1:length(v)
if (isequal(gas,'H2')==1)
    C4=-6.3183e-8;C5=6.33551e-11;
    FJ(k,:) = RJ.*(C1(k)+RJ.*(C2(k)+RJ.*(C3(k)+RJ.*(C4+RJ*C5))));
else
    FJ(k,:) = RJ.*(C1(k)+RJ.*(C2(k)+RJ*C3(k)));
end


GvJ(k,:) = FJ(k,:)+we*VPH(k)-wexe*VPH(k)^2+weye*VPH(k)^3+weze*VPH(k)^4;  %these are the rotational-vibrational term energies in cm-1
end                                                                   %values have been verified vs CARSFT output for N2 and H2 in v = 0.

GvJ1 = circshift(GvJ,-1,1); %shift GvJ array -1 in v index

%%%%%vibrational Raman frequencies in cm-1%%%%%%%%
omegaQ = GvJ1-GvJ;omegaQ = omegaQ(1:size(omegaQ,1)-1,1:length(J)-2);
omegaO = circshift(GvJ1,2,2)-GvJ;omegaO = omegaO(1:size(omegaO,1)-1,3:size(omegaO,2));
omegaS = circshift(GvJ1,-2,2)-GvJ;omegaS = omegaS(1:size(omegaS,1)-1,1:size(omegaS,2)-2);
omegaROT = FJ-circshift(FJ,2,2);omegaROT = omegaROT(1:size(omegaROT,1)-1,3:size(omegaROT,2));


                                                    
                                                 