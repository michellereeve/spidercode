run importLegStrideNormData.m
run importBodyPhaseNormData.m

% pull out chapter 1 intact, flat, ESN data for each species - T/F
bpN_isAranIntact = bpN_isAran & bpN_isIntact;
bpN_isWolfIntactFlatESN = bpN_isWolf & bpN_isIntact & bpN_isFlat & bpN_isESN;
lsN_isAranIntact = lsN_isAran & lsN_isIntact;
lsN_isWolfIntactFlatESN = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN;

% % MIGHT HAVE TO PULL OUT is_trial06 FOR WOLF CHAP1 DATA - CHECK WITH MD

% index back in to pull out actual data (all columns)
bpN_AranIntact = bpN_data(bpN_isAranIntact,:);
bpN_WolfIntactFlatESN = bpN_data(bpN_isWolfIntactFlatESN,:);
lsN_AranIntact = lsN_data(lsN_isAranIntact,:);
lsN_WolfIntactFlatESN = lsN_data(lsN_isWolfIntactFlatESN,:);

% pull out speed data
bpN_AranIntact_speed = cell2mat(bpN_AranIntact(:,9));
bpN_WolfIntactFlatESN_speed = cell2mat(bpN_WolfIntactFlatESN(:,9));
lsN_AranIntact_speed = cell2mat(lsN_AranIntact(:,9));
lsN_WolfIntactFlatESN_speed = cell2mat(lsN_WolfIntactFlatESN(:,9));

% pull out data for bp dataset
bpN_AranIntact_dF = cell2mat(bpN_AranIntact(:,16));
bpN_WolfIntactFlatESN_dF = cell2mat(bpN_WolfIntactFlatESN(:,16)); % dutyFactor
bpN_AranIntact_strideLen = cell2mat(bpN_AranIntact(:,15)); % stride length
bpN_WolfIntactFlatESN_strideLen = cell2mat(bpN_WolfIntactFlatESN(:,15));
bpN_AranIntact_slipF = cell2mat(bpN_AranIntact(:,17)); % stanceSlipFactor
bpN_WolfIntactFlatESN_slipF = cell2mat(bpN_WolfIntactFlatESN(:,17));

% pull out data for ls dataset
lsN_AranIntact_dF = cell2mat(lsN_AranIntact(:,18)); % dutyfactor
lsN_WolfIntactFlatESN_dF = cell2mat(lsN_WolfIntactFlatESN(:,18));

lsN_AranIntact_slipF = cell2mat(lsN_AranIntact(:,23)); % stanceSlipFactor
lsN_WolfIntactFlatESN_slipF = cell2mat(lsN_WolfIntactFlatESN(:,23));

lsN_AranIntact_strideP = cell2mat(lsN_AranIntact(:,15)); % stride Period
lsN_WolfIntactFlatESN_strideP = cell2mat(lsN_WolfIntactFlatESN(:,15));

lsN_AranIntact_stanceP = cell2mat(lsN_AranIntact(:,17)); % stance Period
lsN_WolfIntactFlatESN_stanceP = cell2mat(lsN_WolfIntactFlatESN(:,17));

lsN_AranIntact_swingP = cell2mat(lsN_AranIntact(:,16)); % swing Period
lsN_WolfIntactFlatESN_swingP = cell2mat(lsN_WolfIntactFlatESN(:,16));

lsN_AranIntact_stanceXExc = cell2mat(lsN_AranIntact(:,19)); % stance X excursion (fore-aft)
lsN_WolfIntactFlatESN_stanceXExc = cell2mat(lsN_WolfIntactFlatESN(:,19));

lsN_AranIntact_stanceYExc = cell2mat(lsN_AranIntact(:,21)); % stance Y excursion (med-lat)
lsN_WolfIntactFlatESN_stanceYExc = cell2mat(lsN_WolfIntactFlatESN(:,21));

lsN_AranIntact_swingXExc = cell2mat(lsN_AranIntact(:,20)); % swing X excursion (fore-aft)
lsN_WolfIntactFlatESN_swingXExc = cell2mat(lsN_WolfIntactFlatESN(:,20));

lsN_AranIntact_swingYExc = cell2mat(lsN_AranIntact(:,22)); % swing Y excursion (med-lat)
lsN_WolfIntactFlatESN_swingYExc = cell2mat(lsN_WolfIntactFlatESN(:,22));

% % example polyfit
% [fitresult, gof, output] = createPolyFit(bpN_AranIntact_speed, bpN_AranIntact_dF)
% % click on residuals plot
% title('resid - Aran intact'); % then save & close resid plot
% % click on fit plot
% gtext('Gradient = XX') % copy & paste from output (p1)
% xlabel('StrideAveVel norm');
% ylabel('duty factor');
% title('Polyfit - Aran intact');


%% when I was trying to exclude Aran0022 to figure out what was wrong
% % pull out all Aran Intact except subject 0022
% bpN_isAranIntactNot0022 = cell2mat(bpN_data(:,3))~=22 & bpN_isAran & bpN_isIntact;
% bpN_AranIntactNot0022 = bpN_data(bpN_isAranIntactNot0022,:);
% lsN_isAranIntactNot0022 = cell2mat(lsN_data(:,3))~=22 & lsN_isAran & lsN_isIntact;
% lsN_AranIntactNot0022 = lsN_data(lsN_isAranIntactNot0022,:);

% bpN_AranIntactNot0022_speed = cell2mat(bpN_AranIntactNot0022(:,9));
% lsN_AranIntactNot0022_speed = cell2mat(lsN_AranIntactNot0022(:,9));
% 
% bpN_AranIntactNot0022_dF = cell2mat(bpN_AranIntactNot0022(:,16));
% bpN_AranIntactNot0022_strideLen = cell2mat(bpN_AranIntactNot0022(:,15)); % stride length
% bpN_AranIntactNot0022_slipF = cell2mat(bpN_AranIntactNot0022(:,17)); % stanceSlipFactor
% lsN_AranIntactNot0022_dF = cell2mat(lsN_AranIntactNot0022(:,18)); % dutyfactor
% lsN_AranIntactNot0022_slipF = cell2mat(lsN_AranIntactNot0022(:,23)); % stanceSlipFactor
% lsN_AranIntactNot0022_strideP = cell2mat(lsN_AranIntactNot0022(:,15)); % stride Period
% lsN_AranIntactNot0022_stanceP = cell2mat(lsN_AranIntactNot0022(:,17)); % stance Period
% lsN_AranIntactNot0022_swingP = cell2mat(lsN_AranIntactNot0022(:,16)); % swing Period
% lsN_AranIntactNot0022_stanceXExc = cell2mat(lsN_AranIntactNot0022(:,19)); % stance X excursion (fore-aft)
% lsN_AranIntactNot0022_stanceYExc = cell2mat(lsN_AranIntactNot0022(:,21)); % stance Y excursion (med-lat)
% lsN_AranIntactNot0022_swingXExc = cell2mat(lsN_AranIntactNot0022(:,20)); % swing X excursion (fore-aft)
% lsN_AranIntactNot0022_swingYExc = cell2mat(lsN_AranIntactNot0022(:,22)); % swing Y excursion (med-lat)


