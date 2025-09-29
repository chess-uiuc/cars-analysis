function [xfit,theory,VS,residual] = fitter_interp_3D(fit_save,fit_plot,Tg,Pg)

global data omega_grid
global CARS wlib Plib Tlib DTAUlib

%data=data/abs(trapz(omega_grid,data));
data = data/max(data);

if (size(data,2)==1)
    data=transpose(data);   %this part just makes sure that the 
end                         %data are row vectors that will work 
                            %with the plotting routine
if (size(omega_grid,2)==1)
    omega_grid=transpose(omega_grid);
end

T1=Tg;Tlb=min(Tlib)+1;Tub=max(Tlib);
%P1=Pg;Plb=0.99*P1;Pub=1.01*P1;
P1 = Pg;Plb = min(Plib);Pub = max(Plib);
%if (Plb > Pub);Pa = Plb;Plb = Pub; Pub = Pa;end;clear('Pa');
%if (Plb < min(Plib)); P1 = min(Plib)+1;Plb = 1.0001*P1;Pub = 0.99999*P1;end
ho = 0;

x = [T1 P1 ho 1.0 0.0];   %initial guess
lb = [Tlb Plb 0.9999999*ho 0.999999 -20e-4];  %lower bounds
ub = [Tub Pub 1.0000001*ho 1.000001 20e-4];   %upper bounds

%tic
[xfit,resnorm,residual,exitflag,output,lambda,jacobian]=autogencode_3D_interp(x,lb,ub,2);
%toc

wplot=xfit(4)*(omega_grid-xfit(3));
xfit(6)=resnorm;
P0=1;
theory=interpolate3d_lagrange(CARS,wlib,Tlib,Plib,xfit(1),xfit(2),wplot);
%theory=theory/abs(trapz(wplot,theory));
theory=theory/max(theory);
%theory=theory.*NRB;
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
    annotation('textbox',[0.15, 0.58, 0, 0],'String',['\tau_{23} = ' num2str(xfit(2),'%5.2f')],'FontSize',30,'EdgeColor','none','FitBoxToText','on');
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