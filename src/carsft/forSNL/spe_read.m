function [image] = spe_read(fid,nx,ny,type,point)

%reads in binary image data from PI CCD cameras
%direc = directory path for the image data
%name = name of the .spe file
%nx,ny = pixel format of the .spe file
%type = format and precision of the data
%       LII image accumulations w/ background subtracted are 'float32' or
%       ??
%point = pointer location: = 1 if first image of sequence (header is read)
%                           .ne. 1 otherwise

if point == 1
    header = fread(fid,2050,'int16');
end

image=fread(fid,[nx,ny],type);
image=transpose(image);


        
    
