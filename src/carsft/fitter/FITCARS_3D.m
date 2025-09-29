function [x,resnorm,residual,exitflag,output,lambda,jacobian] = FITCARS_3D(x0,lb,ub,normflag,fitvarnames)

global data omega_grid NRB
global CARS wlib Plib Tlib Alib 
global normfac

if (normflag(1:3) == 'max'); %%normalize the data by max intensity or to unit area
    normfac{1} = max(data);
else
    normfac{1} = abs(trapz(omega_grid,data));
end

normfac{2} = normflag;
data = data/normfac{1};

if (size(data,2)==1)
    data=transpose(data);   %this part just makes sure that the 
end                         %data are row vectors that will work 
                            %with the plotting routine
if (size(omega_grid,2)==1)
    omega_grid=transpose(omega_grid);
end

if (isempty(NRB) == 1);    %if no NRB vector is supplied then 
    NRB = ones(1,length(omega_grid));  %default is to fit NRB-normalized data
else                                    %otherwise multiply theory by NRB
    fprintf('NRB ok.....');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an auto generated MATLAB file from Optimization Tool.
% Start with the default options
PrecondBandWidth_Data = 2;
options = optimset;
% Modify options setting
options = optimset(options,'Display', 'iter-detailed');
options = optimset(options,'Algorithm', 'trust-region-reflective');
options = optimset(options,'MaxIter', 30);
options = optimset(options,'PrecondBandWidth', PrecondBandWidth_Data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@CARS_resid_interp_3D,x0,lb,ub,options); %this part calls the optimization routine%

x(length(x)+1) = resnorm; %put the chi^2 metric in the last column of the fit output vector

fit_plot = 1;fit_save = 0;
plot_fit_results(x,residual,normfac,fitvarnames); %plot the fitresult

end

%%%%%%%%%%%%%%%%%SUBROUTINES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = CARS_resid_interp_3D(x);

%%%%This script computes the residual for a 3-parameter (T, CO, XNR) fit to
%%%%CARS spectra using calls to full frequency-domain calculation
%%%%This version multiplies the theory by the measured argon NRB spectrum
global data omega_grid NRB
global wlib CARS Tlib Plib Alib

T=x(1);
v2=x(2);
v3=x(3);
hshift=x(4);
wexp=x(5);
vshift=x(6);

theory = interpolate3d_lagrange(CARS,wlib,Tlib,Plib,Alib,T,v2,v3,wexp*(omega_grid-hshift));

theory=theory.*NRB;
theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end

indx=find(data ~= 0);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));  %CHANGE
%THIS BACK FOR REAL DATA!!
%f=((data(indx)-vshift)-theory(indx))./sqrt(abs(theory(indx)+1e-5));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_fit_results(xfit,residual,normfac,fitvarnames)

global omega_grid data
global wlib CARS Tlib Plib Alib
global NRB

stg2 = char(fitvarnames{1});
stg3 = char(fitvarnames{2});

xfit

%%%%%%generate the theoretical best-fit result%%%%%
wplot=xfit(5)*(omega_grid-xfit(4));
theory=interpolate3d_lagrange(CARS,wlib,Tlib,Plib,Alib,xfit(1),xfit(2),xfit(3),wplot);
theory=theory.*NRB;

if (normfac{2}(1:3) == 'max')
    theory=theory/max(theory);
    offset = 0.1;
else
    theory=theory/trapz(wplot,theory);
    offset = 0.01;
end

residual=data-(theory+xfit(6));residual=transpose(residual);

%%%%%plot the fit result if required%%%% 
    clf
    hold on
    plot(omega_grid,data,'b','linewidth',4);
    plot(omega_grid,theory+xfit(6),'r--','linewidth',4);
    plot(omega_grid,residual-offset,'-k','linewidth',4);
    plot([min(omega_grid) max(omega_grid)],[-offset -offset],'--k','linewidth',2)
    hold off
    legend('Experiment','Theoretical','Residual','location','NorthWest')
    legend boxoff
    xlabel('Raman Shift (cm^{-1})','Fontsize',30);
    ylabel('CARS Intensity','FontSize',30);
    xlim([min(omega_grid) max(omega_grid)]);
    ylim([1.2*(min(residual)-offset) 1.2*max(data)]);
    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'off'      , ...
      'XColor'      , [0 0 0 ], ...
      'YColor'      , [0 0 0], ...
      'Fontsize'    ,         30, ...
      'LineWidth'   , 3         );
    set(gcf,'color','w');
    
    annotation('textbox',[0.15, 0.65, 0, 0],'String',['T = ' num2str(xfit(1),'%5.2f') ' K'],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    annotation('textbox',[0.15, 0.58, 0, 0],'String',[stg2 ' = ' num2str(xfit(2),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    annotation('textbox',[0.15, 0.51, 0, 0],'String',[stg3 ' = ' num2str(xfit(3),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    %annotation('textbox',[0.15, 0.44, 0, 0],'String',['\chi^{2} = ' num2str(xfit(7),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    
end 