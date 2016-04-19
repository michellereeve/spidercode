%% pull in FILTERED CSVs - need to filter them manually first
run importLegStrideData.m
run importBodyPhaseData.m
%% import arrays of data for subject number + masslength scaling metrics

% imports as individual column vectors, not a cell array
run importScalingMetrics.m

% column vectors:
% Subject
% Lo_in_metres
% NormVel_denom
% NormTime_denom
% NormLength_denom
% NormFreq_denom

%% convert relevant columns by dividing by denominator, for each subject

% which columns correspond to which denominator/column in the
% MassLengthScalingMetrics matrix?
% MassLengthScalingMetrics(:,2) = Lo
% MassLengthScalingMetrics(:,3) = NormVel_denom
% MassLengthScalingMetrics(:,4) = NormTime_denom
% MassLengthScalingMetrics(:,5) = NormFreq_denom
% MassLengthScalingMetrics(:,6) = NormLength_denom

% ls_data(:,9) = strideAveVel = NormVel_denom
% ls_data(:,12) = strideDeltaVel = NormVel_denom
% ls_data(:,15) = stridePeriod = NormTime_denom
% ls_data(:,16) = swingPeriod = NormTime_denom
% ls_data(:,17) = stancePeriod = NormTime_denom
% ls_data(:,19) = stance_x_Excur = NormLength_denom
% ls_data(:,20) = swing_x_Excur = NormLength_denom
% ls_data(:,21) = stance_y_Excur = NormLength_denom
% ls_data(:,22) = swing_x_Excur = NormLength_denom

% bp_data(:,9) = strideAveVel = NormVel_denom
% bp_data(:,10) = strideDeltaVel = NormVel_denom
% bp_data(:,15) = strideLength = NormLength_denom

% for each subject in bp_data and ls_data, columnX =
% columnXdata/relevantDenomValue

%subjectNum = [22 23 26 27 28 29 30 34 36 50 51 52 56 57 58 59 60 61]; %i=18

ls_data_norm = ls_data; % make copies so you can compare old and new
bp_data_norm = bp_data;

for i=1:18
    c_Subject = Subject(i); % current subject
    ls_c_SubjectTF = cell2mat(ls_data(:,3))==c_Subject; % create TF array for current Subject
    bp_c_SubjectTF = cell2mat(bp_data(:,3))==c_Subject; 
    lsdata_c_Subject = ls_data(ls_c_SubjectTF,:);
    bpdata_c_Subject = bp_data(bp_c_SubjectTF,:);
    
    % begin converting things for this subject
    
    % ls_data(:,9) = strideAveVel
    originalData = ls_data(ls_c_SubjectTF,9); % strideAveVel
    newData = (cell2mat(originalData))/NormVel_denom(i); % do the conversion
    newData = num2cell(newData); % convert it from vector into a cell
    ls_data_norm(ls_c_SubjectTF,9) = newData; % replace the correct rows for this subject, in this column, with the converted data
    clear originalData
    clear newData
    
    %ls_data(:,12) = strideDeltaVel
    originalData = ls_data(ls_c_SubjectTF,12);
    newData = (cell2mat(originalData))/NormVel_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,12) = newData;
    clear originalData
    clear newData
    
    %ls_data(:,15) = stridePeriod
    originalData = ls_data(ls_c_SubjectTF,15); 
    newData = (cell2mat(originalData))/NormTime_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,15) = newData;
    clear originalData
    clear newData
    
    % ls_data(:,16) = swingPeriod
    originalData = ls_data(ls_c_SubjectTF,16);
    newData = (cell2mat(originalData))/NormTime_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,16) = newData;
    clear originalData
    clear newData
    
    %ls_data(:,17) = stancePeriod
    originalData = ls_data(ls_c_SubjectTF,17); 
    newData = (cell2mat(originalData))/NormTime_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,17) = newData;
    clear originalData
    clear newData
    
    % ls_data(:,19) = stance_x_Excur
    originalData = ls_data(ls_c_SubjectTF,19); 
    newData = (cell2mat(originalData))/NormLength_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,19) = newData;
    clear originalData
    clear newData
    
    % ls_data(:,20) = swing_x_Excur
    originalData = ls_data(ls_c_SubjectTF,20); 
    newData = (cell2mat(originalData))/NormLength_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,20) = newData;
    clear originalData
    clear newData
    
    % ls_data(:,21) = stance_y_Excur
    originalData = ls_data(ls_c_SubjectTF,21); 
    newData = (cell2mat(originalData))/NormLength_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,21) = newData;
    clear originalData
    clear newData
    
    % ls_data(:,22) = swing_x_Excur
    originalData = ls_data(ls_c_SubjectTF,22);
    newData = (cell2mat(originalData))/NormLength_denom(i);
    newData = num2cell(newData);
    ls_data_norm(ls_c_SubjectTF,22) = newData;
    clear originalData
    clear newData
    
    % bp_data(:,9) = strideAveVel
    originalData = bp_data(bp_c_SubjectTF,9); 
    newData = (cell2mat(originalData))/NormVel_denom(i);
    newData = num2cell(newData);
    bp_data_norm(bp_c_SubjectTF,9) = newData;
    clear originalData
    clear newData
    
    % bp_data(:,10) = strideDeltaVel
    originalData = bp_data(bp_c_SubjectTF,10); 
    newData = (cell2mat(originalData))/NormVel_denom(i);
    newData = num2cell(newData);
    bp_data_norm(bp_c_SubjectTF,10) = newData;
    clear originalData
    clear newData
    
    % bp_data(:,15) = strideLength
    originalData = bp_data(bp_c_SubjectTF,15); 
    newData = (cell2mat(originalData))/NormLength_denom(i);
    newData = num2cell(newData);
    bp_data_norm(bp_c_SubjectTF,15) = newData;
    clear originalData
    clear newData
end


%% save out new CSVs with normalised data

legs_strides_norm_Headers = {'fileName','Species','SubjectNum', 'AblationCond','Terrain','EggsacCond', 'legNumber', 'strideNumber', 'strideAveVel_norm','strideAveVelAng','strideAveYawAng',  ...
    'strideDeltaVel_norm','strideDeltaVelAng','strideDeltaYawAng', 'stridePeriod_norm', 'swingPeriod_norm', ...
    'stancePeriod_norm','dutyFactor','stance_x_Excur_norm','swing_x_Excur_norm','stance_y_Excur_norm','swing_y_Excur_norm','stanceSlipFactor','SegNum','R/L','TrialType'};

% create empty cell array for spreadsheet to save
ls_data_norm2save = cell(size(ls_data_norm,1)+1,size(ls_data_norm,2));

% create empty cell array one row long for headers
ls_data_norm_headersrow = cell(1,size(ls_data_norm,2));

%Create a headers row for the top of the spreadsheet:
[ls_data_norm_headersrow{1,1:end}] = deal(legs_strides_norm_Headers{:});

%Put the headers row together with the data:
ls_data_norm2save = vertcat(ls_data_norm_headersrow,ls_data_norm);

%ls_data_norm = num2cell(ls_data_Norm);

%Save the cell array to a file:
cell2csv('_Legs_Strides_Filtered_Norm.csv', ls_data_norm2save);

body_phase_norm_Headers = {'fileName','Species','SubjectNum', 'AblationCond','Terrain','EggsacCond', 'strideNumber','RefLegFootOffTimes',...
    'strideAveVel_norm', 'strideDeltaVel_norm',...
    'strideAveVelAng','strideDeltaVelAng',...
    'strideAveYawAng', 'strideDeltaYawAng', ...
    'strideLength_norm','dutyFactor', 'stanceSlipFactor',...
    'legPhaseDiffL1','legPhaseDiffL2', 'legPhaseDiffL3','legPhaseDiffL4','legPhaseDiffR1','legPhaseDiffR2','legPhaseDiffR3','legPhaseDiffR4','TrialType'};

% create empty cell array for spreadsheet to save
bp_data_norm2save = cell(size(bp_data_norm,1)+1,size(bp_data_norm,2));

% create empty cell array one row long for headers
bp_data_norm_headersrow = cell(1,size(bp_data_norm,2));

%Create a headers row for the top of the spreadsheet:
[bp_data_norm_headersrow{1,1:end}] = deal(body_phase_norm_Headers{:});

%Put the headers row together with the data:
bp_data_norm2save = vertcat(bp_data_norm_headersrow,bp_data_norm);

%Save the cell array to a file:
cell2csv('_Body_Phase_Filtered_Norm.csv', bp_data_norm2save);