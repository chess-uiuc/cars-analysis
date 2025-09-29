function [xfit,theory,wplot,VS,residual] = FITCARS_2D(fit_save,fit_plot,Tg,Pg,stg2,normflag)

global data omega_grid
global CARS wlib Plib Tlib DTAUlib normfac

if (normflag(1:3) == 'max'); 
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

T1=Tg;Tlb=min(Tlib)+1;Tub=max(Tlib);
P1 = Pg;Plb = min(Plib);Pub = max(Plib);
%P1=Pg;Plb=0.99999*P1;Pub=1.00001*P1;
%P1 = Pg;Plb = 214;Pub = 215;
%if (Plb > Pub);Pa = Plb;Plb = Pub; Pub = Pa;end;clear('Pa');
%if (Plb < min(Plib)); P1 = min(Plib)+1;Plb = 1.0001*P1;Pub = 0.99999*P1;end
ho = 0;

x = [T1 P1 ho 1.0 0.0];   %initial guess
lb = [Tlb Plb -1e-9 0.9999999 0];  %lower bounds
ub = [Tub Pub 1e-9 1.0000001 20e-4];   %upper bounds

%tic
[xfit,resnorm,residual,exitflag,output,lambda,jacobian]=autogencode_2D_interp(x,lb,ub,2);
%toc

wplot=xfit(4)*(omega_grid-xfit(3));
xfit(6)=resnorm;
P0=1;
theory=interpolate2d_lagrange(CARS,wlib,Tlib,Plib,xfit(1),xfit(2),wplot);
if (normfac{2}(1:3) == 'max')
    theory=theory/max(theory);
else
    theory=theory/trapz(wplot,theory);
end

residual=data-(theory+xfit(5));residual=transpose(residual);
%cars_fit_plot(wplot,[data; theory+xfit(5); residual-0.2]);xlim([min(omega_grid) max(omega_grid)]);
pause(0.2);
VS = xfit(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adding in my own plotting features
if fit_plot == 1  
    %figure('Units', 'pixels','Position', [100 100 1000 600]);
    clf
    hold on
        plot(omega_grid,data,'b','linewidth',4);
        plot(omega_grid,theory+xfit(5),'r--','linewidth',4);
        plot(omega_grid,residual-0.01,'k','linewidth',4);
        plot([min(omega_grid) max(omega_grid)],[-0.01 -0.01],'--k','linewidth',2)
    hold off
    legend('Experiment','Theoretical','Residual','location','NorthWest')
    legend boxoff
    xlabel('Raman Shift (cm^{-1})','Fontsize',30);
    ylabel('CARS Intensity','FontSize',30);
    xlim([min(omega_grid) max(omega_grid)]);ylim([1.2*(min(residual)-0.01) 1.2*max(data)]);
    %text(80,0.3,['T = ' num2str(xfit(1)) ' K'],'FontSize',30);
    %text(80,0.3,['\tau = ' num2str(DTAUlib) ' ps'],'FontSize',30);
    
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
    %annotation('textbox',[0.15, 0.58, 0, 0],'String',['\chi^{2} = ' num2str(xfit(6),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    %annotation('textbox',[0.32, 0.68, 0, 0],'String',['\chi^{2} = ' num2str(xfit(6),'%5.3f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    %annotation('textbox',[0.15, 0.58, 0, 0],'String',['\alpha = ' num2str(xfit(2),'%5.2f') ' ps'],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
    
    if fit_save == 1
       set(gcf, 'PaperPositionMode','auto');
       file_place = fullfile('Something');
       print('-dpng','-r300',file_place);  
    end
    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end