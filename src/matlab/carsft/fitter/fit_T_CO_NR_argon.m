function [xfit,resnorm] = fit_T_CO_NR_argon(sqrtflag)

global data omega_grid
global wlib CARS Tlib Plib Alib
global NRB %this version fits theory multiplied by argon nrb spectrum

T1 = 4500;T2 = min(Tlib);T3 = max(Tlib);
C1 = 0.15;C2 = 0; C3 = max(Plib);
alpha1 = 0.05;alpha2 = 0;alpha3 = max(Alib);

h1 = 0; h2 = -10; h3 = 10;
w1 = 1; w2 = 0.95;w3 = 1.05;
v1 = 0; v2 = -0.02*max(data); v3 = -v2;

x0 = [T1 C1 alpha1 h1 w1 v1];
lb = [T2 C2 alpha2 h2 w2 v2];
ub = [T3 C3 alpha3 h3 w3 v3];

data = abs(data/trapz(omega_grid,data));

if (sqrtflag == 'y')
    data = sqrt(abs(data));
end
    
[xfit,resnorm] = fitter_CARS_resid_interp(x0,lb,ub,2,sqrtflag);
xfit = [xfit resnorm];

S = interpolate3d_lagrange(CARS,wlib,Tlib,Plib,Alib,xfit(1),xfit(2),xfit(3),xfit(5)*(omega_grid-xfit(4)));
S = S.*NRB;
S = S/trapz(omega_grid,S);S = abs(S);

if (sqrtflag == 'y');S = sqrt(abs(S));end

S=S-xfit(6);
plot_fit_results(xfit,S,sqrtflag);

end

function plot_fit_results(x,S,sqrtflag)

  global data omega_grid

  P = 1;dwp = 1.4;dtp = 8000;dtau = 0;

  if (size(data,2) ~= size(S,2)); S = transpose(S); end

  residual = S-data-0.1*max(S);
  
  clf;
  plot(omega_grid,data-x(6),omega_grid,S,'--r',omega_grid,residual,'-k','LineWidth',3);
  set(gcf,'units','normalized');
  annotation('textbox',[0.75, 0.90, 0, 0],'String',['T = ' num2str(x(1),4) ' K'],'FontSize',25,'EdgeColor','none','FitBoxToText','on');
  annotation('textbox',[0.75, 0.82, 0, 0],'String',['CO/N_2 = ' num2str(x(2),2)],'FontSize',25,'EdgeColor','none','FitBoxToText','on');
  annotation('textbox',[0.75, 0.74, 0, 0],'String',['\alpha = ' num2str(x(3),2)],'FontSize',25,'EdgeColor','none','FitBoxToText','on');

  set(gca,'LineWidth',2,'XMinorTick','on','YMinorTick','on','FontSize',25);
  xlabel('Raman Shift (cm^{-1})');
  xlim([min(omega_grid)-5 max(omega_grid)+5]);%ylim([1.2*min(residual) 1.05*max(S)]);
  lg = legend('data','theory','residual');
  lg.Position = [0.2 0.8 0 0];
  lg.Box = 'off';
  
  if (sqrtflag == 'y')
      ylabel('(CARS Intensity)^(1/2)');
  else
      ylabel('CARS Intensity');
  end
  
  pause(0.01);

end