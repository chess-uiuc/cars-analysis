%close all
%clear all
clc

%addpath('/Users/spkearn/Documents/Technical/Ultrafast/VCARS_Codes/fscars_v4/vcars_mat/');
% Script to call the Fast CARS program
    % Inputs
    %       omega_grid  = viewing grid in wavenumbers (i.e. 50:0.1:300 cm-1)
    %       T           = temperature of the molecule in question (K)
    %       O2          = O2/N2 number fraction 
    %       dtau_probe  = delay of the probe beam (order of picoseconds)
    %       dwp         = slit width of the spectrometer (i.e. 1 or 1.2 cm-1)
    %       probe_name  = file name of the probe beam temporal profile
    % Outputs
    %       cars        = [1,length(experimental_omega)], CARS signal intensity  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % INPUTS
   omega_grid  = 2100:.05:2380;      
   Tr          = 1000:500:7000;
   Tv          = 1000:500:7000;
   P           = 0.1;  % 0.82 for atm in ABQ
   alpha       = [0:50:1000];
   X           = [1];
   %dtau_probe  = [0];
   %dwp = [dwp fac*dwp];
   dwp = 2.25;
   dtau       = 0;
   dtp         = 8000;
   
    %lib_name = ['Test_library_' num2str(dtau_probe) 'ps.mat'];
    %lib_name = 'impulsive_library_dtp150_dwp1_T_CHINR.mat';
    %lib_name = 'NEQ_Library_04_27_23_voigt_maxnorm.mat';
    %lib_name = 'NEQ_Library_Shot495_1200lmm_voight_maxnorm.mat';
    %lib_name = 'hot_library_for_emccd.mat';
    lib_name = 'PX_NEQ_library_dwp_225_II.mat';
    save_library = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calls the spectrum generating program

%cars = zeros(length(T),length(alpha),length(dwp),length(omega_grid));
cars = zeros(length(Tr),length(Tv),length(alpha),length(omega_grid));
%formatSpec = 'maximum rotational quantum number for %s = %2u\n';
%formatSpec = 'T = %d K, CO = %f, alpha = %f \n';
formatSpec = 'Tr = %d K, Tv = %d,  \alpha %d \n';
%formatSpec = 'T = %d K, alpha = %d, dwp = %d cm-1 \n';

for k = 1:length(Tr)
    for j = 1:length(Tv)
        for l = 1:length(alpha)
       fprintf(formatSpec,Tr(k), Tv(j), alpha(l));
       %cars(k,j,:) = picosecond_CARS_v4a(omega_grid,T(k),P,dtau,dtaus,2,dts,dtp,dwp,alpha(j));
       cars(k,j,l,:) = CARSFT_dev(omega_grid,[Tr(k) Tv(j)],P,X,dtp,dtau,alpha(l),dwp);
       cars(k,j,l,:) = cars(k,j,l,:)/max(cars(k,j,l,:));
       %cars(k,j,l,:) = CARSFT_dev(omega_grid,T(k),P,X,dtp,dtau,alpha(j),dwp(l));
       %cars(k,j,l,:) = cars(k,j,l,:)/max(cars(k,j,l,:));
       plot(omega_grid,squeeze(cars(k,j,l,:)),'LineWidth',3);
       ylim([-0.05 1.05]);xlim([min(omega_grid) max(omega_grid)]);
       set(gca,'FontSize',25);set(gca,'LineWidth',3);
       set(gca,'XMinorTick','on','YMinorTick','on');
       xlabel('Raman Shift (cm^{-1})');ylabel('CARS Intensity');
       pause(0.001);
       end
    end
end


    % Save the library
    if save_library == 1
        % dwp, dtau, Tlib, CARS,wlib
            Tlib = Tr;
            Plib = Tv;
            Alib = alpha;
            wlib = omega_grid;
            CARS = cars;
        
        save(['computed_libraries/' lib_name],'CARS','Tlib','Plib','Alib','wlib','dwp','-mat');
        %save(['computed_libraries/' lib_name],'CARS','Tlib','Plib','wlib','dwp','P','dts','dtp','dtau','-mat');
    end
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% 
% if spectral_plot == 1
% 
% % Plotting 
% cc = hsv(5);
% figure('Units', 'pixels','Position', [100 100 775 575]);
%     hold on
%         plot(omega_grid,cars/max(cars),'Color',cc(1,:),'linewidth',2);
%     hold off
%     xlabel('Raman Shift (cm$^{-1}$)','interpreter','Latex','Fontsize',10);
%     ylabel('CARS Intensity','interpreter','Latex');
%     text(190,0.7,['P = ' num2str(P) ' atm'],'FontSize',20);
%     text(190,0.6,['T = ' num2str(T) ' K'],'FontSize',20);
%     %text(260,0.9,['O$_2$/N$_2$ = ' num2str(O2)],'interpreter','Latex','FontSize',15);
%     text(190,0.5,['$\tau$ = ' num2str(dtau_probe) ' ps'],'interpreter','Latex','FontSize',20);
%     text(190,0.4,['O$_2$/N$_2$ = ' num2str(O2)],'interpreter','Latex','FontSize',20);
%     xlim([0 300])
%     set(gca, ...
%       'Box'         , 'off'     , ...
%       'TickDir'     , 'out'     , ...
%       'TickLength'  , [.02 .02] , ...
%       'XMinorTick'  , 'on'      , ...
%       'YMinorTick'  , 'on'      , ...
%       'YGrid'       , 'off'      , ...
%       'XColor'      , [0 0 0 ], ...
%       'YColor'      , [0 0 0], ...
%       'Fontsize'    ,         25, ...
%       'LineWidth'   , 3         );
%     set(gcf,'color','w');
%     %set(gca,'xtick',[],'ytick',[]);
%     %axis off
%     
% end
   