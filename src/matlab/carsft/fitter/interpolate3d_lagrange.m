function [G] = interpolate3d_lagrange(f,w,x,y,z,xo,yo,zo,wo)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function takes a library of 1-D functions, f, that are
% parameterized by X and Y (library grid locations)
% and returns an function, fo, evaluated at the query points
% (xo,yo) via bilinear interpolation
% the result is further interpolated onto the coarse grid, wo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Find the closest indices for the 3D library to the input guesses
[g,indx]=sort(abs(xo-x));indx=indx(1:4);indx=sort(indx);  
[g,indy]=sort(abs(yo-y));indy=indy(1:4);indy=sort(indy);
[g,indz]=sort(abs(zo-z));indz=indz(1:4);indz=sort(indz);

T(1) = x(indx(1)); T(2) = x(indx(2)); T(3) = x(indx(3)); T(4) = x(indx(4));
O(1) = y(indy(1)); O(2) = y(indy(2)); O(3) = y(indy(3)); O(4) = y(indy(4));
d(1) = z(indz(1)); d(2) = z(indz(2)); d(3) = z(indz(3)); d(4) = z(indz(4));


    % Generate the weight terms for lagrange interpolation
    T1 = ((xo-T(2))*(xo-T(3))*(xo-T(4)))/((T(1)-T(2))*(T(1)-T(3))*(T(1)-T(4)));
    T2 = ((xo-T(1))*(xo-T(3))*(xo-T(4)))/((T(2)-T(1))*(T(2)-T(3))*(T(2)-T(4)));
    T3 = ((xo-T(1))*(xo-T(2))*(xo-T(4)))/((T(3)-T(1))*(T(3)-T(2))*(T(3)-T(4)));
    T4 = ((xo-T(1))*(xo-T(2))*(xo-T(3)))/((T(4)-T(1))*(T(4)-T(2))*(T(4)-T(3)));
    O1 = ((yo-O(2))*(yo-O(3))*(yo-O(4)))/((O(1)-O(2))*(O(1)-O(3))*(O(1)-O(4)));
    O2 = ((yo-O(1))*(yo-O(3))*(yo-O(4)))/((O(2)-O(1))*(O(2)-O(3))*(O(2)-O(4)));
    O3 = ((yo-O(1))*(yo-O(2))*(yo-O(4)))/((O(3)-O(1))*(O(3)-O(2))*(O(3)-O(4)));
    O4 = ((yo-O(1))*(yo-O(2))*(yo-O(3)))/((O(4)-O(1))*(O(4)-O(2))*(O(4)-O(3)));
    d1 = ((zo-d(2))*(zo-d(3))*(zo-d(4)))/((d(1)-d(2))*(d(1)-d(3))*(d(1)-d(4)));
    d2 = ((zo-d(1))*(zo-d(3))*(zo-d(4)))/((d(2)-d(1))*(d(2)-d(3))*(d(2)-d(4)));
    d3 = ((zo-d(1))*(zo-d(2))*(zo-d(4)))/((d(3)-d(1))*(d(3)-d(2))*(d(3)-d(4)));
    d4 = ((zo-d(1))*(zo-d(2))*(zo-d(3)))/((d(4)-d(1))*(d(4)-d(2))*(d(4)-d(3)));


% confirm that we are 'in the box'
 %   [x(indx(1)) xo x(indx(2))]
  %  [y(indy(1)) yo y(indy(2))]
  %  [z(indz(1)) zo z(indz(2))]
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate at the FIRST slit width
    f111(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(1),indz(1),:))),wo); 
    f121(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(2),indz(1),:))),wo); 
    f131(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(3),indz(1),:))),wo); 
    f141(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(4),indz(1),:))),wo); 
    f211(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(1),indz(1),:))),wo);
    f221(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(2),indz(1),:))),wo);
    f231(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(3),indz(1),:))),wo);
    f241(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(4),indz(1),:))),wo);
    f311(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(1),indz(1),:))),wo);
    f321(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(2),indz(1),:))),wo);
    f331(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(3),indz(1),:))),wo);
    f341(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(4),indz(1),:))),wo);
    f411(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(1),indz(1),:))),wo);
    f421(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(2),indz(1),:))),wo);
    f431(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(3),indz(1),:))),wo);
    f441(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(4),indz(1),:))),wo);
    
    F=[f111 f121 f131 f141;f211 f221 f231 f241; f311 f321 f331 f341;f411 f421 f431 f441];

    g1=zeros(1,1,size(F,3)); %pre-allocate memory for the output vector

    for k=1:size(F,3)
        g1(1,1,k)=[T1 T2 T3 T4]*F(:,:,k)*[O1;O2;O3;O4];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Interpolate at the SECOND slit width
    f112(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(1),indz(2),:))),wo); 
    f122(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(2),indz(2),:))),wo); 
    f132(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(3),indz(2),:))),wo); 
    f142(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(4),indz(2),:))),wo); 
    f212(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(1),indz(2),:))),wo);
    f222(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(2),indz(2),:))),wo);
    f232(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(3),indz(2),:))),wo);
    f242(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(4),indz(2),:))),wo);
    f312(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(1),indz(2),:))),wo);
    f322(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(2),indz(2),:))),wo);
    f332(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(3),indz(2),:))),wo);
    f342(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(4),indz(2),:))),wo);
    f412(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(1),indz(2),:))),wo);
    f422(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(2),indz(2),:))),wo);
    f432(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(3),indz(2),:))),wo);
    f442(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(4),indz(2),:))),wo);
    
    F=[f112 f122 f132 f142;f212 f222 f232 f242; f312 f322 f332 f342;f412 f422 f432 f442];
    


    g2=zeros(1,1,size(F,3)); %pre-allocate memory for the output vector

    for k=1:size(F,3)
        g2(1,1,k)=[T1 T2 T3 T4]*F(:,:,k)*[O1;O2;O3;O4];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate at the THIRD slit width
    f113(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(1),indz(3),:))),wo); 
    f123(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(2),indz(3),:))),wo); 
    f133(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(3),indz(3),:))),wo); 
    f143(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(4),indz(3),:))),wo); 
    f213(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(1),indz(3),:))),wo);
    f223(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(2),indz(3),:))),wo);
    f233(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(3),indz(3),:))),wo);
    f243(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(4),indz(3),:))),wo);
    f313(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(1),indz(3),:))),wo);
    f323(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(2),indz(3),:))),wo);
    f333(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(3),indz(3),:))),wo);
    f343(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(4),indz(3),:))),wo);
    f413(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(1),indz(3),:))),wo);
    f423(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(2),indz(3),:))),wo);
    f433(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(3),indz(3),:))),wo);
    f443(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(4),indz(3),:))),wo);
    
    F=[f113 f123 f133 f143;f213 f223 f233 f243; f313 f323 f333 f343;f413 f423 f433 f443];
    

    g3=zeros(1,1,size(F,3)); %pre-allocate memory for the output vector

    for k=1:size(F,3)
        g3(1,1,k)=[T1 T2 T3 T4]*F(:,:,k)*[O1;O2;O3;O4];
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate at the FOURTH slit width
    f114(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(1),indz(4),:))),wo); 
    f124(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(2),indz(4),:))),wo); 
    f134(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(3),indz(4),:))),wo); 
    f144(1,1,:)=interp1(w,squeeze(squeeze(f(indx(1),indy(4),indz(4),:))),wo); 
    f214(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(1),indz(4),:))),wo);
    f224(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(2),indz(4),:))),wo);
    f234(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(3),indz(4),:))),wo);
    f244(1,1,:)=interp1(w,squeeze(squeeze(f(indx(2),indy(4),indz(4),:))),wo);
    f314(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(1),indz(4),:))),wo);
    f324(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(2),indz(4),:))),wo);
    f334(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(3),indz(4),:))),wo);
    f344(1,1,:)=interp1(w,squeeze(squeeze(f(indx(3),indy(4),indz(4),:))),wo);
    f414(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(1),indz(4),:))),wo);
    f424(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(2),indz(4),:))),wo);
    f434(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(3),indz(4),:))),wo);
    f444(1,1,:)=interp1(w,squeeze(squeeze(f(indx(4),indy(4),indz(4),:))),wo);
    
    F=[f114 f124 f134 f144;f214 f224 f234 f244; f314 f324 f334 f344;f414 f424 f434 f444];
    

    g4=zeros(1,1,size(F,3)); %pre-allocate memory for the output vector

    for k=1:size(F,3)
        g4(1,1,k)=[T1 T2 T3 T4]*F(:,:,k)*[O1;O2;O3;O4];
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Interpolate the two results between the slit width to get the actual resulting theory spectrum

    F = [g1;g2;g3;g4];
    
    
    G=zeros(1,size(F,3));
    
    for k = 1:size(F,3)
        G(k) = [d1 d2 d3 d4]*F(:,:,k);
    end

%G=G/(dx*dy*dz);
%fo=interp1(w,fo,wo);
G=G/max(G);
