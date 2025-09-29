function F = rot_terms(J,v)

%this function calculates rotational term values of N2 in cm-1 

%%%%molecular constants FOR N2 pulled directly from CARS.MOL%%%%%
Be = [0.19982600E+01 0.14456220E+01];
alphae = [0.17303500E-01 0.15932680E-01] ;
game = [-0.31536099E-04 0.000000];
De = [0.57740000E-05 0.57740000E-05];
betae = [0.15500000E-07 0.000000];
delte = [0.0000000 0.000000];
ho = [0.30000000E-11 0.000000];
he = [0.18000000E-11 0.000000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




