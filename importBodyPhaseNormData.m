%% Import data from text file - BodyPhase filtered normalised data
% Script for importing data from the following text file:
%
%    /Users/Michelle/Research/data/__DATA_ANALYSIS/_ALL_DATA/_Results/_Body_Phase_Filtered_Norm.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/04/18 20:34:35

%% Initialize variables.
filename = '/Users/Michelle/Research/data/__DATA_ANALYSIS/_ALL_DATA/_Results/_Body_Phase_Filtered_Norm.csv';
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
BodyPhaseFilteredNorm = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

bpN_data = BodyPhaseFilteredNorm;
clear BodyPhaseFilteredNorm  

bpN_isWolf = strcmp(bpN_data(:,2),'Wolf');
bpN_isAran = strcmp(bpN_data(:,2),'Aran');
bpN_isFlat = strcmp(bpN_data(:,5),'F');
bpN_isRough = strcmp(bpN_data(:,5),'R');
bpN_isESN = strcmp(bpN_data(:,6),'ESN');
bpN_isESY = strcmp(bpN_data(:,6),'ESY');
bpN_isIntact = strcmp(bpN_data(:,4),'intact');
bpN_isR4ablt = strcmp(bpN_data(:,4),'R4ablt');
bpN_isL3aR4m = strcmp(bpN_data(:,4),'L3aR4m');