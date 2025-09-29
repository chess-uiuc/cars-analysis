function [G] = interpolate2d_lagrange(f,w,x,y,xo,yo,wo)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function takes a library of 1-D functions, f, that are
% parameterized by X and Y (library grid locations)
% and returns an function, fo, evaluated at the query points
% (xo,yo) via largange multipliers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Find the closest indices for the 3D library to the input guesses
[g,indx]=sort(abs(xo-x));indx=indx(1:4);indx=sort(indx);  
[g,indy]=sort(abs(yo-y));indy=indy(1:4);indy=sort(indy);

T(1) = x(indx(1)); T(2) = x(indx(2)); T(3) = x(indx(3)); T(4) = x(indx(4));
O(1) = y(indy(1)); O(2) = y(indy(2)); O(3) = y(indy(3)); O(4) = y(indy(4));


    % Generate the weight terms for lagrange interpolation
    T1 = ((xo-T(2))*(xo-T(3))*(xo-T(4)))/((T(1)-T(2))*(T(1)-T(3))*(T(1)-T(4)));
    T2 = ((xo-T(1))*(xo-T(3))*(xo-T(4)))/((T(2)-T(1))*(T(2)-T(3))*(T(2)-T(4)));
    T3 = ((xo-T(1))*(xo-T(2))*(xo-T(4)))/((T(3)-T(1))*(T(3)-T(2))*(T(3)-T(4)));
    T4 = ((xo-T(1))*(xo-T(2))*(xo-T(3)))/((T(4)-T(1))*(T(4)-T(2))*(T(4)-T(3)));
    O1 = ((yo-O(2))*(yo-O(3))*(yo-O(4)))/((O(1)-O(2))*(O(1)-O(3))*(O(1)-O(4)));
    O2 = ((yo-O(1))*(yo-O(3))*(yo-O(4)))/((O(2)-O(1))*(O(2)-O(3))*(O(2)-O(4)));
    O3 = ((yo-O(1))*(yo-O(2))*(yo-O(4)))/((O(3)-O(1))*(O(3)-O(2))*(O(3)-O(4)));
    O4 = ((yo-O(1))*(yo-O(2))*(yo-O(3)))/((O(4)-O(1))*(O(4)-O(2))*(O(4)-O(3)));

   O1;
   O2;
   O3;
   O4;

    
% confirm that we are 'in the box'
   % [x(indx(1))  x(indx(2)) xo x(indx(3)) x(indx(4))]
    %[y(indy(1)) yo y(indy(2))]
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate at the FIRST slit width
    f111(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(1),:))),wo); 
    f121(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(2),:))),wo); 
    f131(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(3),:))),wo); 
    f141(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(4),:))),wo); 
    f211(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(1),:))),wo);
    f221(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(2),:))),wo);
    f231(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(3),:))),wo);
    f241(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(4),:))),wo);
    f311(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(1),:))),wo);
    f321(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(2),:))),wo);
    f331(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(3),:))),wo);
    f341(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(4),:))),wo);
    f411(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(1),:))),wo);
    f421(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(2),:))),wo);
    f431(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(3),:))),wo);
    f441(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(4),:))),wo);
    
    F=[f111 f121 f131 f141;f211 f221 f231 f241; f311 f321 f331 f341;f411 f421 f431 f441];

    g1=zeros(1,1,size(F,3)); %pre-allocate memory for the output vector

    for k=1:size(F,3)
        g1(1,1,k)=[T1 T2 T3 T4]*F(:,:,k)*[O1;O2;O3;O4];
    end
%size(F)

    G = squeeze(g1); G = G';

%G=G/(dx*dy*dz);
%fo=interp1(w,fo,wo);
G=G/max(G);
