%% LEG STRIDE FILTERED DATA FILE
%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/Michelle/Research/data/__DATA_ANALYSIS/_ALL_DATA/_Results/LegStridesCompiledDataFiltered.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/03/23 15:34:27

% Initialize variables.
filename = '/Users/Michelle/Research/data/__DATA_ANALYSIS/_ALL_DATA/_Results/LegStridesCompiledDataFiltered.csv';
delimiter = ',';
startRow = 2;

% Format string for each line of text:
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
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
dataArray([3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]) = cellfun(@(x) num2cell(x), dataArray([3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), 'UniformOutput', false);
LegStridesCompiledDataFiltered = [dataArray{1:end-1}];
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% SELECT SUBSETS OF DATA

ls_data = LegStridesCompiledDataFiltered;
clear LegStridesCompiledDataFiltered

% % select subsets of data - this way is more risky as 'ismember' takes more
% % inputs but might be harder to figure out if it's working properly for
% % that reasons
% ls_WolfInd = find(ismember(ls_data(:,2),'Wolf'));
% ls_AranInd = find(ismember(ls_data(:,2),'Aran'));
% ls_WolfIntactInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact'));
% ls_AranIntactInd = find(ismember(ls_data(:,2),'Aran') & ismember(ls_data(:,4),'intact'));
% ls_WolfR4abltInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'R4ablt'));
% ls_AranR4abltInd = find(ismember(ls_data(:,2),'Aran') & ismember(ls_data(:,4),'R4ablt'));
% ls_WolfL3aR4mInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'L3aR4m'));
% ls_AranL3aR4mInd = find(ismember(ls_data(:,2),'Aran') & ismember(ls_data(:,4),'L3aR4m'));
% ls_AllIntact = find(ismember(ls_data(:,4),'intact'));
% ls_AllR4ablt = find(ismember(ls_data(:,4),'R4ablt'));
% ls_AllL3aR4m = find(ismember(ls_data(:,4),'L3aR4m'));
% ls_WolfIntactFlatInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact') & ismember(ls_data(:,5),'F'));
% ls_WolfIntactFlatESNInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact') & ismember(ls_data(:,5),'F') & ismember(ls_data(:,6),'ESN'));
% ls_WolfIntactFlatESYInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact') & ismember(ls_data(:,5),'F') & ismember(ls_data(:,6),'ESY'));
% ls_WolfIntactRougInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact') & ismember(ls_data(:,5),'R'));
% ls_WolfIntactRougESNInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact') & ismember(ls_data(:,5),'R') & ismember(ls_data(:,6),'ESN'));
% ls_WolfIntactRougESYInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'intact') & ismember(ls_data(:,5),'R') & ismember(ls_data(:,6),'ESY'));
% ls_WolfR4abltFlatInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'R4ablt') & ismember(ls_data(:,5),'F'));
% ls_WolfR4abltRougInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'R4ablt') & ismember(ls_data(:,5),'R'));
% ls_WolfL3aR4mFlatInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'L3aR4m') & ismember(ls_data(:,5),'F'));
% ls_WolfL3aR4mRougInd = find(ismember(ls_data(:,2),'Wolf') & ismember(ls_data(:,4),'L3aR4m') & ismember(ls_data(:,5),'R'));

% select for subsets - gives 0 or 1 False/True (these are all cols of
% strings)
ls_isWolf = strcmp(ls_data(:,2),'Wolf');
ls_isAran = strcmp(ls_data(:,2),'Aran');
ls_isFlat = strcmp(ls_data(:,5),'F');
ls_isRough = strcmp(ls_data(:,5),'R');
ls_isESN = strcmp(ls_data(:,6),'ESN');
ls_isESY = strcmp(ls_data(:,6),'ESY');
ls_isIntact = strcmp(ls_data(:,4),'intact');
ls_isR4ablt = strcmp(ls_data(:,4),'R4ablt');
ls_isL3aR4m = strcmp(ls_data(:,4),'L3aR4m');

% need to treat columns of numbers differently
ls_isSeg0 = cell2mat(ls_data(:,24))==0;
ls_isSeg1 = cell2mat(ls_data(:,24))==1;
ls_isSeg2 =  cell2mat(ls_data(:,24))==2;
ls_isSeg3 =  cell2mat(ls_data(:,24))==3;

% you can combine these subsets like so:
ls_isAranFlatIntact_SegNum0 = ls_isAran & ls_isFlat & ls_isIntact & ls_isSeg0;

% and then index back into the original cell array of data to pull out just
% the subset you want
AranIntactFlat_Seg0 = ls_data(ls_isAranFlatIntact_SegNum0,:);


