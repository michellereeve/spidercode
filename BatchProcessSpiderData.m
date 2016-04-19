
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

%Create a new folder to save the results files into
mkdir([pathName '_Results' sep]);
saveDir = dirCheck([pathName '_Results' sep ]);

strides_CompiledCellArray = cell(1,26);
bodyPhase_CompiledCellArray = cell(1,26);

for f = 1:length(fileNamesArray)
    
c_fileName = fileNamesArray{f};

[bp_compiledCellArray,body_phase_Headers, ls_compiledCellArray, legs_strides_Header, filePrefix, combinedEvents] = ProcessSpiderData(c_fileName,pathName);
close all;

if f == 1
    [strides_CompiledCellArray{1,1:end}] = deal(legs_strides_Header{:});
    [bodyPhase_CompiledCellArray{1,1:end}] = deal(body_phase_Headers{:});
end

strides_CompiledCellArray = vertcat(strides_CompiledCellArray,ls_compiledCellArray); %#ok<AGROW>
bodyPhase_CompiledCellArray = vertcat(bodyPhase_CompiledCellArray,bp_compiledCellArray); %#ok<AGROW>

end

%Create a new 'Output' folder to save the processed data files into


 cell2csv([saveDir '_LegsStridesCompiledData.csv' ], strides_CompiledCellArray);
 cell2csv([saveDir '_BodyPhaseCompiledData.csv' ], bodyPhase_CompiledCellArray);
 
 

 
