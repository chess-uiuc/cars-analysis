function [x,resnorm] = FITCARS(x0,lb,ub,N,normflag,fitvarnames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a more generalized fitting code using interpolated spectral libraries
%
% INPUTS
% x0 = initial guess vector
% lb,ub = lower and upper bounds
% N = number of library parameters (2 or 3 for now)
% normflag = 'max' to normalize data to unit max, any other string for area
%               normalized data
% fitvarnames = cell array of strings containing legend strings for library
%           parameters 2 (2D library) or 2 and 3 (3D library). Default 
%           label for parameter 1 is Temperature in Kelvin units
%
% OUTPUTS
% x = vector of best-fit results with chi^2 (resnorm) appended in last
%       column
% resnorm = chi^2 vector, as specified in user-supplied objective function
%              to lsqnonlin routine in Matlab optimization toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%this part calls the optimization routine
if N == 2;
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
       lsqnonlin(@CARS_resid_interp_2D,x0,lb,ub,options);
elseif N == 3;
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
       lsqnonlin(@CARS_resid_interp_3D,x0,lb,ub,options); %this part calls the optimization routine%
else
    fprintf('number of library parameters is not consistent')
    return
end

x(length(x)+1) = resnorm; %put the chi^2 metric in the last column of the fit output vector

fit_plot = 1;fit_save = 0;
plot_fit_results(x,residual,normfac,fitvarnames,N); %plot the fitresult

end

%%%%%%%%%%%%%%%%%SUBROUTINES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = CARS_resid_interp_3D(x);

%%%%This script computes the residual for a 3-parameter (e.g. T, CO, XNR) fit to
%%%%CARS spectra using calls to the LaGrange interpolated 3D library
%%%%This version multiplies the theory by the measured argon NRB spectrum
%%%%Simply make NRB = ones(1,length(omega_grid)) to fit to argon-normalized
%%%%data

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
function f = CARS_resid_interp_2D(x);
%%%%This script computes the residual for a 2-parameter (e.g. T, XNR) fit to
%%%%CARS spectra using calls to the LaGrange interpolated 2D library
%%%%This version multiplies the theory by the measured argon NRB spectrum
%%%%Simply make NRB = ones(1,length(omega_grid)) to fit to argon-normalized
%%%%data

global data omega_grid NRB
global wlib CARS Tlib Plib Alib

T=x(1);
v2=x(2);
hshift=x(3);
wexp=x(4);
vshift=x(5);

theory = interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,v2,wexp*(omega_grid-hshift));
theory=theory.*NRB;
theory=theory/abs(trapz(wexp*(omega_grid-hshift),theory));

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end
indx=find(data ~= 0);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_fit_results(xfit,residual,normfac,fitvarnames,N)

global omega_grid data
global wlib CARS Tlib Plib Alib
global NRB


%%%%%%generate the theoretical best-fit result%%%%%
if N == 2
    stg2 = char(fitvarnames);
    HS = xfit(3);WEXP = xfit(4);VS = xfit(5);
    wplot=WEXP*(omega_grid-HS);
    theory=interpolate2d_lagrange(CARS,wlib,Tlib,Plib,xfit(1),xfit(2),wplot);
end

if N == 3
    stg2 = char(fitvarnames{1}); stg3 = char(fitvarnames{2});
    HS = xfit(4);WEXP = xfit(5);VS = xfit(6);
    wplot=WEXP*(omega_grid-HS);
    theory=interpolate3d_lagrange(CARS,wlib,Tlib,Plib,Alib,xfit(1),xfit(2),xfit(3),wplot);
end

theory=theory.*NRB;

if (normfac{2}(1:3) == 'max')
    theory=theory/max(theory);
    offset = 0.1;
else
    theory=theory/abs(trapz(wplot,theory));
    offset = 0.01;
end

residual=data-(theory+VS);residual=transpose(residual);
offset = max(residual)+0.01;

%%%%%plot the fit result if required%%%% 
    clf
    hold on
    plot(omega_grid,data,'b','linewidth',4);
    plot(omega_grid,theory+VS,'r--','linewidth',4);
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
    
    if N == 3
        annotation('textbox',[0.15, 0.51, 0, 0],'String',[stg3 ' = ' num2str(xfit(3),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    end
    
    %annotation('textbox',[0.15, 0.44, 0, 0],'String',['\chi^{2} = ' num2str(xfit(7),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    
end 