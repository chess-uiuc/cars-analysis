function [fo] = interpolate2d(f,w,x,y,xo,yo,wo)

% INPUTS
% f = 1-D function values of library 
% w = 1-D coordinate (x-axis or wavenumber) array of library
% x = fitting variable 1 (usually temperature) array of library
% y = fitting variable 2 array of library
% xo = guess for variable 1
% yo = guess for variable 2
% wo = wavenumber grid * expansion value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function takes a library of 1-D functions, f, that are
% parameterized by X and Y (library grid locations)
% and returns an function, fo, evaluated at the query points
% (xo,yo) via bilinear interpolation
% the result is further interpolated onto the coarse grid, wo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[g,indx]=sort(abs(xo-x));indx=indx(1:2);indx=sort(indx);  %find the (i,j) indices of the nearest two (x,y) locations to (xo,yo)
[g,indy]=sort(abs(yo-y));indy=indy(1:2);indy=sort(indy);

% confirm that we are 'in the box'
    %[x(indx(1)) xo x(indx(2))]
    %[y(indy(1)) yo y(indy(2))]

dx=x(indx(2))-x(indx(1));dy=y(indy(2))-y(indy(1));  %calculate the grid size (dx,dy)

%f11=f(indx(1),indy(1),:);f12=f(indx(1),indy(2),:);  %extract the library values for all z at the 4 closest (x,y) points
%f22=f(indx(2),indy(2),:);f21=f(indx(2),indy(1),:);

f11(1,:)=squeeze(f(indx(1),indy(1),:)); 
f12(1,:)=squeeze(f(indx(1),indy(2),:));
f22(1,:)=squeeze(f(indx(2),indy(2),:));
f21(1,:)=squeeze(f(indx(2),indy(1),:));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


g11(1,1,:)=squeeze(interp1(w,f11,wo)); %extract the library values for all z at the 4 closest (x,y) points
g12(1,1,:)=squeeze(interp1(w,f12,wo)); %and interpolate them onto the experimental grid
g22(1,1,:)=squeeze(interp1(w,f22,wo));% ??should we renormalize these to unit maximum??
g21(1,1,:)=squeeze(interp1(w,f21,wo));


F=[g11 g12;g21 g22]; %form Nz pages of 2x2 matrices containing the values of f at the 4 closest points

 DX=[x(indx(2))-xo xo-x(indx(1))];
 DY=[y(indy(2))-yo;yo-y(indy(1))];  %form the weights for linear interpolation

fo=zeros(1,size(F,3)); %pre-allocate memory for the output vector


for k=1:size(F,3)
    fo(k)=[x(indx(2))-xo xo-x(indx(1))]*F(:,:,k)*[y(indy(2))-yo;yo-y(indy(1))];
end

cc = f11;
%fo=fo/(dx*dy);

%fo=interp1(w,fo,wo);

fo=fo/max(fo); max(fo);
