
%Batch process spider data and save to CSV file

%Utility variable
if isempty(strfind(computer,'MAC'))
    sep = '\';
else sep = '/';
end

%Get information about file names,  ask  user browse to folder
dataFolder = uigetdir(cd, 'Browse to recording data folder with digitised files');
cd(dataFolder)
pathName = [dataFolder sep];
dirStructure = dir([pathName '*.csv']);
[fileNamesArray{1:size(dirStructure)}] = deal(dirStructure.name);

strides_CompiledCellArray = cell(1,12);
bodyPhase_CompiledCellArray = cell(1,16);

for f = 1:length(fileNamesArray)
    
c_fileName = fileNamesArray{f};

%[bp_compiledCellArray,body_phase_Headers, ls_compiledCellArray, legs_strides_Header, filePrefix, combinedEvents] = ProcessSpiderData_Chloe(c_fileName,pathName);
[bp_compiledCellArray,body_phase_Headers, ls_compiledCellArray, legs_strides_Header, filePrefix, combinedEvents] = ProcessSpiderData_OrbW(c_fileName,pathName);
close all;

if f == 1
    [strides_CompiledCellArray{1,1:end}] = deal(legs_strides_Header{:});
    [bodyPhase_CompiledCellArray{1,1:end}] = deal(body_phase_Headers{:});
end

strides_CompiledCellArray = vertcat(strides_CompiledCellArray,ls_compiledCellArray); %#ok<AGROW>
bodyPhase_CompiledCellArray = vertcat(bodyPhase_CompiledCellArray,bp_compiledCellArray); %#ok<AGROW>


end

 cell2csv([pathName 'LegStridesCompiledData.csv' ], strides_CompiledCellArray)
 cell2csv([pathName 'BodyPhaseCompiledData.csv' ], bodyPhase_CompiledCellArray)

 