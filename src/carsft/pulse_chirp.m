t =  [-10000:.1:10000];
dt = 40;
tau = 20;
f = exp(-2*log(2)*((t-tau)/dt).^2);

N = 2;
dw = N*14.7/dt;
DT = dt*1000;
beta = (dw/(2*dt)).^2 - ((1e12*log(2))/(pi*3e10*dt*dt)).^2;
beta = sqrt(beta);
%beta = beta*1000;

co = 3e10;
phi = beta*2*pi*co*1e-12*(t-tau).*(t-tau);

g = f.*exp(-i*phi);

F = fft(f);G = fft(g);
F = fftshift(F);G = fftshift(G);
w = time_to_freq(t*1000,2)/3e10;

F = F/sqrt(max(F.*conj(F)));G = G/sqrt(max(G.*conj(G)));
plot(w,F.*conj(F),w,G.*conj(G));xlim([-2 2]);

ff = ifft(F); gg = ifft(G);
ff = ff/sqrt(max(ff.*conj(ff)));gg = gg/sqrt(max(gg.*conj(gg)));
plot(t,ff.*conj(ff),t,gg.*conj(gg));
t2 = freq_to_time(w*3e10)/1000;
