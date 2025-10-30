kb = 0.6950356; %Boltmann constant in cm-1/K
T = 1000; %tempertature (K)

f = read_O2_props(fname,15);
g = f.g; F = f.term; v = f.v; J = f.J; N = f.N;


Nb = g.*exp(-F/(kb*T));
nssw = mod(N,2); %spin statistical weight (=1 for Nodd, = 0 for Neven)
Nb = nssw.*Nb; Nb = Nb/sum(Nb);


%isotropic Raman Q branch frequencies: DN = 0, DJ = 0;

