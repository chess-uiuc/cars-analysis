function Tv = area_Tv(wlims,omega_grid,cars,NRB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to compute Tv from the ratio of v = 0 and v = 1
% CARS signals
%
% wlims = limits of integration (2 X 2) row 1 is v = 0 and row 2 is v = 1
% omega_grid = experimental wavenumber grid
% CARS = experimental CARS spectrum (NOT NORMALIZED by the NRB)
% NRB = argon NRB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
we = 2358.518; %we in cm-1
theta = 1.44*we; %characteristic vibrartional temperature (K)

% find pixel indices for v = 0,1,and 2 peaks
v0 = find(omega_grid >= wlims(1,1) & omega_grid <= wlims(1,2));
v1 = find(omega_grid >= wlims(2,1) & omega_grid <= wlims(2,2));
v2 = find(omega_grid >= wlims(3,1) & omega_grid <= wlims(3,2));

nback = 5;
back0 = [ [min(v0)-nback:1:min(v0)-1] [max(v0)+1:1:max(v0)+nback] ]; %take local backround at 
back1 = [ [min(v1)-nback:1:min(v1)-1] [max(v1)+1:1:max(v1)+nback] ]; %+/- nback pixels about wlims
back2 = [ [min(v2)-nback:1:min(v2)-1] [max(v2)+1:1:max(v2)+2] ];

if (size(wlims,1) == 4) %add the v = 3 peak if needed must have at least v = 0-2 to perform this analysis
    v3 = find(omega_grid >= wlims(4,1) & omega_grid <= wlims(4,2));
    back3 = [ [min(v3)-nback:1:min(v3)-1] [max(v3)+1:1:max(v3)+2] ];
end

cars = bsxfun(@rdivide,cars,NRB); %normalize the data by the NRB reference spectrum

n = length(omega_grid);
Tv = zeros(size(cars,1)-1,1);
for k = 1:length(Tv)
    c = mean(cars(k:k+1,:)); %running average fit
    back = mean(c(n-20:n)); %estimate background from last 20 pixels in spectrum
    c = c - back;

    a0 = c(v0)-mean(c(back0));a0 = abs(trapz(omega_grid(v0),a0));a0=sqrt(a0);
    a1 = c(v1)-mean(c(back1));a1 = abs(trapz(omega_grid(v1),a1));a1=sqrt(a1)/2; 
    a2 = c(v2)-mean(c(back2));a2 = abs(trapz(omega_grid(v2),a2));a2=sqrt(a2)/3;

    if (size(wlims,1)==3)
        f = (a0+a1+a2)/(a1+a2);
    else
        a3 = c(v3)-mean(c(back3));a3 = abs(trapz(omega_grid(v3),a3));a3=sqrt(a3)/4;
        f = (a0+a1+a2+a3)/(a1+a2+a3);
    end
    Tv(k) = theta/log(f);

    clf;
    plot(omega_grid,c,'LineWidth',1);hold on
    plot(omega_grid,zeros(1,length(omega_grid)),'--k');
    plot(omega_grid(v0),c(v0)-mean(c(back0)),'-r','Linewidth',3);
    plot(omega_grid(v1),c(v1)-mean(c(back1)),'-r','LineWidth',3);
    plot(omega_grid(v2),c(v2)-mean(c(back2)),'-r','LineWidth',3);

    if (size(wlims,1) == 4)
        plot(omega_grid(v3),c(v3)-mean(c(back3)),'-r','LineWidth',3);
    end

    set(gca,'XMinorTick','on','YMinorTick','on');set(gca,'FontSize',25);set(gca,'LineWidth',3);
    xlabel('Raman Shift (cm^{-1})');
    an1 = annotation('textbox',[0.2 0.6 0 0],'String',['T_{v} = ' num2str(Tv(k),4) ' K'],...
        'FontSize',20,'FitBoxToText','on','LineStyle','none');
    hold off
    pause(0.2)

end





