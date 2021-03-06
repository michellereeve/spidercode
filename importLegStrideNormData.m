%% Import data from text file. LegsStrides, filtered normalised data
% Script for importing data from the following text file:
%
%    /Users/Michelle/Research/data/__DATA_ANALYSIS/_ALL_DATA/_Results/_Legs_Strides_Filtered_Norm.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/04/18 20:33:16

%% Initialize variables.
filename = '/Users/Michelle/Research/data/__DATA_ANALYSIS/_ALL_DATA/_Results/_Legs_Strides_Filtered_Norm.csv';
delimiter = ',';
startRow = 2;

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: double (%f)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
dataArray([3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]) = cellfun(@(x) num2cell(x), dataArray([3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), 'UniformOutput', false);
LegsStridesFilteredNorm = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

lsN_data = LegsStridesFilteredNorm;
clear LegsStridesFilteredNorm

lsN_isWolf = strcmp(lsN_data(:,2),'Wolf');
lsN_isAran = strcmp(lsN_data(:,2),'Aran');
lsN_isFlat = strcmp(lsN_data(:,5),'F');
lsN_isRough = strcmp(lsN_data(:,5),'R');
lsN_isESN = strcmp(lsN_data(:,6),'ESN');
lsN_isESY = strcmp(lsN_data(:,6),'ESY');
lsN_isIntact = strcmp(lsN_data(:,4),'intact');
lsN_isR4ablt = strcmp(lsN_data(:,4),'R4ablt');
lsN_isL3aR4m = strcmp(lsN_data(:,4),'L3aR4m');

% need to treat columns of numbers differently
lsN_isSeg0 = cell2mat(lsN_data(:,24))==0;
lsN_isSeg1 = cell2mat(lsN_data(:,24))==1;
lsN_isSeg2 =  cell2mat(lsN_data(:,24))==2;
lsN_isSeg3 =  cell2mat(lsN_data(:,24))==3;

lsN_isLegL1 = cell2mat(lsN_data(:,7))==1;
lsN_isLegL2 = cell2mat(lsN_data(:,7))==2;
lsN_isLegL3 = cell2mat(lsN_data(:,7))==3;
lsN_isLegL4 = cell2mat(lsN_data(:,7))==4;
lsN_isLegR1 = cell2mat(lsN_data(:,7))==5;
lsN_isLegR2 = cell2mat(lsN_data(:,7))==6;
lsN_isLegR3 = cell2mat(lsN_data(:,7))==7;
lsN_isLegR4 = cell2mat(lsN_data(:,7))==8;

% % you can combine these subsets like so:
% lsN_isAranFlatIntact_SegNum0 = lsN_isAran & lsN_isFlat & lsN_isIntact & lsN_isSeg0;
% 
% % and then index back into the original cell array of data to pull out just
% % the subset you want
% AranIntactFlat_Seg0 = lsN_data(lsN_isAranFlatIntact_SegNum0,:);