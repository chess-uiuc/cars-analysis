function N = v_fraction(v,we,T);

%%computes vibrational Boltzmann fraction for a diatomic molecule

hck=1.44; %hc/k in K/cm-1;
v=[0:1:v];

N=zeros(length(v),length(we));
for k=1:length(we)
    N(:,k)=exp(-v*hck*we(k)/T)*(1-exp(-hck*we(k)/T));
end

%N=bsxfun(@rdivide,N,sum(N,1));   %normalizes total population to 1 for each species for small amounts in higher v levels