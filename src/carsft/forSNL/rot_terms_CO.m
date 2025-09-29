function F = rot_terms_CO(J,v)

%this function calculates rotational term values of CO in cm-1 

%read molecular data from CARS.MOL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
molpar = read_cars_dot_mol('CARS.MOL');
%[gas,WE,WX,WY,WZ,ALPHAE,BE,DE,BETAE,DELTE,GAME,H0,HE,gI] = read_cars_dot_mol('CARS.MOL');

kgas = find(molpar.gasname == 'CO');

Be = molpar.BE(kgas);
alphae = molpar.ALPHAE(kgas);
game = molpar.GAME(kgas);
De = molpar.DE(kgas);
betae = molpar.BETAE(kgas);
delte = molpar.DELTE(kgas);
ho = molpar.H0(kgas);
he = molpar.HE(kgas);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RJJ = J.*(J+1);
F=zeros(length(J),length(v),length(Be));

for j=1:length(Be)
for k=1:length(v)
    vph = v(k)+0.5;
    C1 = Be(j) + vph*(-alphae(j) + vph*game(j));
    C2 = De(j) + vph*(betae(j) + vph*delte(j));C2=-C2;
    C3 = ho(j) + vph*he(j);
    
    F(:,k,j) = RJJ.*(C1 + RJJ.*(C2+RJJ*C3));   %rotational term energies in cm-1
end
end




