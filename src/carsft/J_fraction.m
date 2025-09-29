function [NJ] = J_fraction(J,Be,gI,T)

%computes rotational Boltzmann fraction given
%J quantum numbers
%Be in cm-1
%nuclear spin degeneracies for odd and even J, gI
%the temperature, T

hck=1.441304;   %hc/k in K/cm-1
JJ = J; %store input J vector here

%calculate the partition function via direct sum
maxJ = 201;J=[0:1:maxJ]; %assumes all thermal population will be accounted for with Jmax = 200
if (mod(maxJ,2) == 0);maxJ=maxJ+1;end

Jeven = [0:2:maxJ-1];keven = [1:2:maxJ];
Jodd = [1:2:maxJ];kodd = [2:2:maxJ+1];

NJ = zeros(length(J),1);

NJ(keven) = (2*Jeven+1).*exp(-hck*Be*Jeven.*(Jeven+1)/T);
NJ(kodd) = (2*Jodd+1).*exp(-hck*Be*Jodd.*(Jodd+1)/T);
Q = sum(NJ);
NJ = NJ/Q;

%now account for the gI
NJ(keven) = gI(2)*NJ(keven);
NJ(kodd) = gI(1)*NJ(kodd);
NJ = 2*NJ/(sum(gI));

%truncate back to input J values
NJ = NJ(1:length(JJ));
