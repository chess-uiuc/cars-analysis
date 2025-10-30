function t = freq_to_time(f);

%function computes the time array for an inverse DFT from frequency f
%(in sec-1) to time IN FEMTOSECONDS

dt=1/max(abs(2*f));
N=length(f);

t = [0:1:N-1]-floor(N/2);
t = t*dt;

C=1e15;  %conversion factor to convert to fs. You can easily include this as an input and output time in any units
t=t*C;