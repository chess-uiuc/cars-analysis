function [xfit,resnorm] = extract_slit(raman_shift,c,argon,normflag)

global omega_grid data
global cars_params
global flag

%direc = ['/Users/spkearn/Documents/Technical/Hypersonics/Shock_Tunnel_CARS/Run470/'];
%load([direc '2Torr_air_spectrum.mat'],'raman_shift','c','argon');

% c = average cars spectrum at known condition
% raman_shift = wavenumber vector (cm-1)
% argon = average argon spectrum
% flag = 'max' normalize data to unit max, unit area otherwise

flag = normflag;

wmin = 2310;wmax = 2340;
indx = find(raman_shift >= wmin & raman_shift <= wmax);
omega_grid = raman_shift(indx);
data = c; data = data - mean(data(1:100));data = data(indx);

if (isequal(flag,'max'))
    data = data/max(data);
else
    data = data/abs(trapz(omega_grid,data));
end

cars_params = [292 3/760 1 8000 0 0];
%dwp = [1 6.5]; 
%hs = -1.2;wexp = 1.0;
%x = [dwp hs wexp];

%xo = [0.8 0.2 -.5 1];
xo = [1.5 .4 0.25 1];
lb = [0.1 0.1 -4 0.95];
ub = [3 10 4 1.05];

PrecondBandWidth_Data = 2;
options = optimset;
options = optimset(options,'Display', 'iter-detailed');
options = optimset(options,'Algorithm', 'trust-region-reflective');
options = optimset(options,'MaxIter', 30);
options = optimset(options,'PrecondBandWidth', PrecondBandWidth_Data);

[xfit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
       lsqnonlin(@slit_resid,xo,lb,ub,options);

w = xfit(4)*omega_grid-xfit(3);
S = CARSFT_dev(w,cars_params(1),cars_params(2),cars_params(3),cars_params(4),...
    cars_params(5),cars_params(6),[xfit(1) xfit(2)]);

if (isequal(flag,'max'))
    S = S/max(S);
else
    S = S/abs(trapz(w,S));
end

clf; plot(omega_grid,data,omega_grid,S,'--r','LineWidth',3);
set(gca,'XMinorTick','on','YMinorTick','on');set(gca,'FontSize',25);set(gca,'LineWidth',3);
ylabel('CARS Intensity');xlabel('Raman Shift (cm^{-1})');
a1 = annotation('textbox',[0.2, 0.8 0 0],'String',['\sigma = ' num2str(xfit(1),3) ' cm^{-1}'],'FontSize',25);
a1.FitBoxToText = 'on';a1.LineStyle = 'none';
a2 = annotation('textbox',[0.2, 0.7 0 0],'String',['\Gamma = ' num2str(xfit(2),3) ' cm^{-1}'],'FontSize',25);
a2.FitBoxToText = 'on';a2.LineStyle = 'none';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = slit_resid(x)

%this is the objective function to be minimized
%
%INPUTS


global omega_grid data
global cars_params
global flag

%1. parse the input cell array, x
dwp = [x(1) x(2)];
HS = x(3);WEXP = x(4);

cp = cars_params;


%3. calculate the theoretical spectrum
w = WEXP*omega_grid-HS; %shifts the wave# axis for the theory calculation
S = CARSFT_dev(w,cp(1),cp(2),cp(3),cp(4),cp(5),cp(6),dwp);

if (isequal(flag,'max'))
    S = S/max(S);
else
    S = S/abs(trapz(w,S));
end

%4. compute the objective function
f = (data - S)./sqrt(abs(data)+1e-5);

plot(omega_grid,data,omega_grid,S,omega_grid,f-0.1);

end
