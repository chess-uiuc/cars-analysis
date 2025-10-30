function [cfit,waven] = validate_carsft_mat

%point to next folder up to find the current CARSFT_mat version
carsft_direc = pwd;indx = strfind(carsft_direc,'/');indx = indx(length(indx));
carsft_direc = carsft_direc(1:indx);addpath(carsft_direc);
addpath('CARSFT_OG_FEB_2024/');

%%%validation case #1%%%%
%%%nanosecond, T = 1504 K, N2 = 1, CO = 0.3, dwp = 1%%%%%

cfitdataname = 'T1504N250CO15XNRZEROISO.dat';
cfit = import_carsft(cfitdataname);
waven = cfit(:,1);cfit = cfit(:,2);
plot(waven,cfit);

T = 1504;X = [1 0.3];P = 1;
dtp = 8000;dtau = 0;alpha = 0;
dwp = 1;
S = CARSFT_dev(waven,T,P,X,dtp,dtau,alpha,dwp);

clf;
S = sqrt(S);S = S/max(S);
cfit = cfit/max(cfit);
plot(waven,S,waven,cfit,'--r','LineWidth',2);
set(gca,'LineWidth',2,'XMinorTick','on','YMinorTick','on');
set(gca,'FontSize',25);
xlabel('Raman Shift (cm^{-1})');ylabel('(CARS Intensity)^{1/2}');
ylim([-0.02 1.02]);
lg = legend('CRF CARSFT','Matlab');lg.Box = 'off';
lg.Position = [0.15   0.7459    0.1833    0.1571];



function carsft_res = import_carsft(filename, startRow, endRow)

%% Initialize variables.
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%18f%14f%14f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

carsft_res = [dataArray{1:end-1}];
