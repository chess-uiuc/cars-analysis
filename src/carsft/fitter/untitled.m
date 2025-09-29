global omega_grid data
global cars_params

direc = ['/Users/spkearn/Documents/Technical/Hypersonics/Shock_Tunnel_CARS/Run468/'];
load([direc '2Torr_air_spectrum.mat'],'raman_shift','c','argon');

wmin = 2280;wmax = 2360;
indx = find(raman_shift >= wmin & raman_shift <= wmax);
omega_grid = raman_shift(indx);
data = mean(c);data = data - mean(data(1:100));data = data(indx);

cars_params = [292 2/760 1 8000 0 0];
dwp = [1 6.5]; 
hs = -1.2;wexp = 1.0;
x = [dwp hs wexp];

xo = [1 5 -1.2 1];
lb = [0.1 0.1 -4 0.95];
ub = [3 10 4 1.05];

PrecondBandWidth_Data = 2;
options = optimset;
options = optimset(options,'Display', 'iter-detailed');
options = optimset(options,'Algorithm', 'trust-region-reflective');
options = optimset(options,'MaxIter', 30);
options = optimset(options,'PrecondBandWidth', PrecondBandWidth_Data);

[xfit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
       lsqnonlin(@test_function,xo,lb,ub,options);

w = xfit(4)*omega_grid-xfit(3);
S = CARSFT_dev(w,cars_params(1),cars_params(2),cars_params(3),cars_params(4),...
    cars_params(5),cars_params(6),dwp);
S = S/abs(trapz(w,S));
plot(omega_grid,data,omega_grid,S,'--r','LineWidth',3);
set(gca,'XMinorTick','on','YMinorTick','on');set(gca,'FontSize',25);set(gca,'LineWidth',3);
ylabel('CARS Intensity');xlabel('Raman Shift (cm^{-1}');
