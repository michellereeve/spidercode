% Transform leg phases into spreadsheet for circular ANOVA

run PullOutDataForC1Stats.m
clear lsN_is* bpN_is* lsN_Aran* bpN_Aran* c_bpN_Aran* c_lsN_Aran* *FESNC_* *FESNC

saveDir = ['/Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/'];

% bpN_WolfIntactFESNUnC
% what columns do we need?
% columns 18-25 - leg phase diffs
% col 7 - stride number
% col 9 - stride ave vel
% Then need from ls_data (or manual)
% LvsR L=0 R=1
% segNum 0,1,2,3

legPhases_vert = vertcat(c_bpN_WolfIntactFESNUnC_legPhaseDiffL1,c_bpN_WolfIntactFESNUnC_legPhaseDiffL2,c_bpN_WolfIntactFESNUnC_legPhaseDiffL3,...
    c_bpN_WolfIntactFESNUnC_legPhaseDiffL4,c_bpN_WolfIntactFESNUnC_legPhaseDiffR1,c_bpN_WolfIntactFESNUnC_legPhaseDiffR2,c_bpN_WolfIntactFESNUnC_legPhaseDiffR3,...
    c_bpN_WolfIntactFESNUnC_legPhaseDiffR4); 

% the original array bpN_WolfIntactFESNUnC is 51 rows long. So vertically stacking up
% all the legs gives 51x8 = 408 rows total. Each block of 51 is for one
% leg, and as the bp_data array only measured things relative to leg L2,
% all the information can be duplicated. Therefore strideNum & vel are
% going to be the same 51 values repeated 8 times.

vel_vert = repmat(bpN_WolfIntactFESNUnC_vel,8,1); % vertically stack all the strideAveVel - replicate array by 8 times, in the vetical direction

strideNum = cell2mat(bpN_WolfIntactFESNUnC(:,7)); % pull out stride number column
strideNum_vert = repmat(strideNum,8,1); % vert stack all stride nums, replicate 8 times in vert direction

% legNum, LvsR and segNum aren't in this array, but they will follow the
% same pattern - it will be the same 8 values each repeated 51 times e.g.
% legNum1 x51, legNum2 x51... legNum8 x51. All the legs are in a fixed
% order and each block of 51 = a single leg, so segNum and LvsR can be
% duplicated for each element within the block of 51 rows.

legNum = [1;2;3;4;5;6;7;8];
legNum_vert = repelem(legNum,51); % repeat each element in vector 51 times 

LvsR = [0;0;0;0;1;1;1;1];
LvsR_vert = repelem(LvsR,51);

segNum = [0;1;2;3;0;1;2;3];
segNum_vert = repelem(segNum,51);

legPhaseSpreadsheet = horzcat(legNum_vert,legPhases_vert,vel_vert,LvsR_vert,segNum_vert,strideNum_vert); % put together everything
headers={'leg','legPhaseDiff','StrideAveVel','LvsR','segNum','strideNum'};

legPhaseSpreadsheet_h = vertcat(headers,num2cell(legPhaseSpreadsheet));

cell2csv([saveDir,'legphases_forCircANOVA_WolfIntactFESNUnC.csv'], legPhaseSpreadsheet_h);

clear i headers 




