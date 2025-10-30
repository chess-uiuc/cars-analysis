function f = SlitObjFunction(x)

global omega_grid data
global cars_params

% cars_params(1) = temperature (K)
% cars_params(2) = Pressure (atm)
% cars_params(3) = X mole fraction
% cars_params(4) = probe pulse duration, dtp (ps)
% cars_params(5) = probe pulse delay, dtau (ps)
% cars_params(6) = nonresonant alpha factor

dwp = [x(1) x(2)];
hshift = x(3);
wexp = x(4);
vshift = x(5);

w = wexp*omega_grid - hshift;

S = CARSFT_dev(w,cars_params(1),cars_params(2),cars_params(3),cars_params(4),cars_params(5)...
    cars_params(6),dwp);

f = ((data-vshift)-S)./sqrt(abs(data)+1e-6);
