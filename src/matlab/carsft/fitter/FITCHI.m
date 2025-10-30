function [x,resnorm] = FITCHI(x0,lb,ub,N,normflag,fitvarnames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a more generalized fitting code using interpolated spectral libraries
% NOTE: this code requires a library of the complex-valued susceptibility
%        the result is singly comvolved with a voight slit width
%
%  This will work directly for ns CARS with a transform-limited probe,
%  where the probe width is << than the Raman linewidth and there is zero
%  probe delay
% 
%  for fs CARS convolve the probe with phase factor with complex chi and
%  place that result in the library 
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
% GLOBALS
% slit = slit function. This will be interpolated onto the fine wave# grid
%        used in the susceptibility library.  Supply as a 2D array as
%        [islit, wslit].
% CARS = spectral intensity library (NTlib x NPlib x NAlib x Nwlib) for N =
%        3 or (NTlib x NPlib x Nwlib) for N = 2;
% Tlib = First library physics variable. ALmost always rotational
%         temperature (K)
% Plib = Second library physics variable
% Alib = Third library physics variable, for N = 3
% wlib = library wavenumber grid (cm-1)
% omega_grid = experimental wave# array (cm-1)
% data = experimental data. Code assumes NO argon normalization and
%         multiplies theory by the NRB vector
% NRB = NRB argon spectrum. must be same length as data and omega_grid. Can
%        leave empty for argon-normalized data or set equal to vector of
%        ones
% weights = must be same length as data and omega_grid. Vector of weight
%           factors. Eech wave# point is multiplied by weights(k) in the
%           least-squares objective function. Leave empty or set = vector
%           of ones for no preferential weighting.
% normface = 'max' for unit max normalization of data. Anything else for
%               unit area normalization
%
% OUTPUTS
% x = vector of best-fit results with chi^2 (resnorm) appended in last
%       column
% resnorm = chi^2 vector, as specified in user-supplied objective function
%              to lsqnonlin routine in Matlab optimization toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global data omega_grid NRB
global CARS wlib Plib Tlib Alib 
global normfac weights slit

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
end                                    %otherwise multiply theory by NRB

if (isempty(weights) == 1)
    weights = ones(1,length(omega_grid));
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

%this part preparaes the slit function
if(find(size(slit)==2) == 1);slit = transpose(slit);end %work with column data
wcol = find(min(slit) < 0); %this index is the wave# column (has negative entry)
if (wcol == 1);tmp=slit(:,1);slit(:,1)=slit(:,2);slit(:,2)=tmp;clear('tmp');end
wslit = max(abs(slit(:,2)));wslit = [-wslit:mean(diff(wlib)):wslit];
slit = interp1(slit(:,2),slit(:,1),wslit); %interpolate user supplied slit onto susceptibility grid

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
global weights

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
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));
f = weights.*f; %apply optional weights to emphasize some parts of spectrum more

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
global weights slit

T=x(1);
v2=x(2);
hshift=x(3);
wexp=x(4);
vshift=x(5);


theory = interpolate2d_lagrange(CARS,wlib,Tlib,Plib,T,v2,wlib);
theory = conv_spk(slit,theory.*conj(theory));  %NOTE: slit fxn is computed once on unshifted/unstretched grid. this is ok for small wexp 
theory = interp1(wlib,theory,wexp*(omega_grid-hshift));
theory=theory.*NRB;
theory=theory/abs(trapz(omega_grid,theory));
%theory(find(isnan(theory)==1))=0;

if (size(data,1) ~= size(theory,1));theory=transpose(theory);end
indx=find(data ~= 0);
f=((data(indx)-vshift)-theory(indx))./sqrt(abs(data(indx)+1e-5));
f = weights.*f;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = real_voigt(x,sigma,gamma)

%function to compute the real-valued voigt function on 1D grid, x
%
% x = wave# grid
% sigma = gaussian width parameter (std. dev)
% gamma = lorentz width parameter (FWHM)

xv = x/sigma;
yv = gamma/(2*sigma);
z = xv+i*yv;
F = faddeeva(z,16);

f = real(F);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f2 = conv_spk(f,g) %simple covolution routine

%check the length of the inputs to determine the output vector size
%this should be the length of the largest input vector
if (length(f) < length(g)); M = length(g);else; M = length(f);end

%check length of input vector to treat N = even case by adding a zero to
%the end of the vector
if (mod(length(f),2) == 0); f = [f 0];end
if (mod(length(g),2) == 0); g = [g 0];end

N = abs(length(f)-length(g));
if length(f) < length(g)
    f = [zeros(1,N/2) f zeros(1,N/2)];
else
    g = [zeros(1,N/2) g zeros(1,N/2)];
end

F = fft(f);G = fft(g);
f2 = F.*G;f2 = ifftshift(ifft(f2));

f2 = f2(1:M);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_fit_results(xfit,residual,normfac,fitvarnames,N)

global omega_grid data
global wlib CARS Tlib Plib Alib
global NRB slit


%%%%%%generate the theoretical best-fit result%%%%%
if N == 2
    stg2 = char(fitvarnames);
    HS = xfit(3);WEXP = xfit(4);VS = xfit(5);
    wplot=WEXP*(omega_grid-HS);
    theory=interpolate2d_lagrange(CARS,wlib,Tlib,Plib,xfit(1),xfit(2),wlib);
    theory=conv_spk(slit,theory.*conj(theory));
    theory = interp1(wlib,theory,WEXP*(omega_grid-HS));
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
offset = max(residual);

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