%% this script is for testing the effects of cooling on wolf spider speed

% compare trial number 3 (intact, flat, no ES, part of the ES data
% collection -not cooled) and trial number 5 (intact, flat, no ES, part of
% the ablation data collection - cooled)

run PullOutDataForC1Stats.m

% for testing effects of cooling:

% just pull out the subjects that overlap both cooled and non-cooled conditions
bpN_isWolfCoolComp  = cell2mat(bpN_data(:,3))==56 | cell2mat(bpN_data(:,3))==57 | cell2mat(bpN_data(:,3))==58 | cell2mat(bpN_data(:,3))==60 | cell2mat(bpN_data(:,3))==61; 
bpN_isWolfIntactFlatESNCoolComp = bpN_isWolf & bpN_isIntact & bpN_isFlat & bpN_isESN & bpN_isWolfCoolComp;

bpN_WolfIntactFlatESNCoolComp = bpN_data(bpN_isWolfIntactFlatESNCoolComp,:);


% first need speed for wolf intact flat ESN trialtype=3 and trialtype=5

bpN_isWolfIntactFlatESNCoolCompT3 = cell2mat(bpN_WolfIntactFlatESNCoolComp(:,26))==3;
bpN_isWolfIntactFlatESNCoolCompT5 = cell2mat(bpN_WolfIntactFlatESNCoolComp(:,26))==5;

bpN_WolfIntactFlatESNCoolCompT3 = bpN_WolfIntactFlatESNCoolComp(bpN_isWolfIntactFlatESNCoolCompT3,:);
bpN_WolfIntactFlatESNCoolCompT5 = bpN_WolfIntactFlatESNCoolComp(bpN_isWolfIntactFlatESNCoolCompT5,:);

% pull out speed
bpN_WolfIntactFlatESNCoolCompT3_speed = cell2mat(bpN_WolfIntactFlatESNCoolCompT3(:,9));
bpN_WolfIntactFlatESNCoolCompT5_speed = cell2mat(bpN_WolfIntactFlatESNCoolCompT5(:,9));

% run unpaired t-test
[bp_h,bp_p,bp_ci,bp_stats] = ttest2(bpN_WolfIntactFlatESNCoolCompT3_speed,bpN_WolfIntactFlatESNCoolCompT5_speed);

% run Cohen's D for effect sizes
[M1,M2,S1,S2,n1,n2,Stdev_p,Cohens_d] = CohensD(bpN_WolfIntactFlatESNCoolCompT3_speed,bpN_WolfIntactFlatESNCoolCompT5_speed);

% boxplot -- THIS BIT NEEDS TESTING WITHIN THIS SCRIPT
plotData = cell2mat(bpN_WolfIntactFlatESNCoolComp(:,9));
trialGroup = cell2mat(bpN_WolfIntactFlatESNCoolComp(:,26));
boxplot(plotData,trialGroup);
ylabel('StrideAveVel');
xlabel('3 = uncooled, 5 = cooled');

