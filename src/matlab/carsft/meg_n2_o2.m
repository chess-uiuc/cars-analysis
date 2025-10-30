function [GAMMA,gamma,DE] = meg_n2_o2(Ji,T,branch);

To=295;   %reference temperature (K)

%%%%%fit parameters from Farrow et al., Appl. Opt. 26 (1987) in comments%%%%%
%%%%From Lavorel et al. Opt. Lett. 20 (1995 in use%%%%%%%%%%%%%%
alpha=[0.02646 0.0167];%0.0231; %cm-1 atm-1
beta=[1.850 1.45];%1.67;
delta=[1.199 1.32];%1.21;
m=0.1381;%0.1487;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0.69503476;  %Boltzmann's constant in cm-1/K

J=[0:1:2*length(Ji)];  %add additional J values to ensure that sum of off-diagonal matrix elements converges
%J=Ji;

GAMMA=zeros(length(alpha),length(J));  %initialize the output linewidth vector

v=0;  %consider the ground vibrational state

gamma=zeros(length(J),length(J),length(alpha));  %G-matrix initialized page 1 is N2, page 2 is O2 G=G(final,initial)
DE=zeros(length(J));   %change in energy for transition from rotational state j (cols) to state i (rows)

F=rot_terms(J,0);F=squeeze(F);  %rotational term values in cm-1 in v = 0 for the specified J

f(1)=(1-exp(-m))/(1-exp(-m*T/To))*sqrt(To/T);  %pre-exponential factor developed by Farrow et al. (1987)
f(2)=(To/T)^1.32;

for m=1:length(alpha)

for j = 1:length(J)
    DE(:,j)=F(j,m)-F(:,m);   %fill the matrix DEij. The change in term energy (cm-1) for a transtion from state i to state j
end

C1=alpha(m).*f(m);
C1=C1*((1+1.5*F(:,m)/(k*T*delta(m)))./(1+1.5*F(:,m)/(k*T))).^2;

D=exp(-beta(m)*DE/(k*T));

for j=1:length(J)
    for i=1:j-1
        gamma(j,i,m)=C1(i)*D(i,j);   %calculate the matrix elements for energy-increasing collisions (below the diagonal)
    end
end

for j=1:length(J)
    for i=1:j-1
        gamma(i,j,m)=(2*J(i)+1)/(2*J(j)+1)*gamma(j,i)*exp(DE(i,j)/(k*T));  %calculate the matrix elements for energy-decreasing collisions (above the diagonal)
    end
end

ind = find(mod(J,2)==0);  %find J = even collisions
GAMMA(m,ind)=2*sum(gamma(ind,ind,m),1); %factor of 2 for FWHM linewidths

ind = find(mod(J,2)==1); %find J = odd collisions
GAMMA(m,ind)=2*sum(gamma(ind,ind,m),1);  %factor of 2 for FWHM linewidths

%%These are Q-branch linewidths. Make random-phase approximation.
%%This takes the average of Q(J) and Q(J+2) for S-branch linewidths
%%and Q(J) and Q(J-2) for O-branch lines. 

%%OUTPUT S-BRANCH LINEWIDTHS HERE AND USE j-2 SUBSCRIPT IN MAIN CODE. BE
%%CAREFUL!!

%GAMMA=transpose(GAMMA);
%ind = find(mod(J,2)==0);GAMMA(ind,m)=0.5*(GAMMA(ind,m)+circshift(GAMMA(ind,m),-1));
%ind = find(mod(J,2)==1);GAMMA(ind,m)=0.5*(GAMMA(ind,m)+circshift(GAMMA(ind,m),-1));

%GAMMA=transpose(GAMMA);

end

GAMMA=transpose(GAMMA);

GAMMA = GAMMA(1:length(Ji),:);  %tuncate back to input J vector only

gamma=gamma(1:length(Ji),1:length(Ji),:); %reduce gamma back to length of initial J vector for output


%%%%%this part sets values of G matrix for transitions between odd and even J to zero
rownums=[1:1:size(gamma,1)];rownums=transpose(rownums);rownums=repmat(rownums,1,length(rownums));
colnums=transpose(rownums);
R=ones(size(gamma));

[ind1,ind2]=ind2sub(size(gamma),find(mod(colnums,2)==0 & mod(rownums,2)==1));
R(ind1,ind2)=0;

[ind1,ind2]=ind2sub(size(gamma),find(mod(rownums,2)==0 & mod(colnums,2)==1));
R(ind1,ind2)=0;

gamma=R.*gamma;  %set parity odd/even and even/odd values to zero