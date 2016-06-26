% Create rose plots

run PullOutDataForC1Stats.m
clear lsN_is* bpN_is* lsN_Aran* bpN_Aran* c_bpN_Aran* c_lsN_Aran* *FESNC_* *FESNC % clear all non- uncooled stuff and all aran

[bpN_WolfIntactFESNUnC_SLOW, bpN_WolfIntactFESNUnC_MED, bpN_WolfIntactFESNUnC_FAST] = SetSpeedBins(bpN_WolfIntactFESNUnC,bpN_WolfIntactFESNUnC_vel);

% calculate median values for overlay onto rose plots
legPhaseDiffs_rads_s = circ_ang2rad(cell2mat(bpN_WolfIntactFESNUnC_SLOW(:,18:25)));
legPhaseDiffs_rads_m = circ_ang2rad(cell2mat(bpN_WolfIntactFESNUnC_MED(:,18:25)));
legPhaseDiffs_rads_f = circ_ang2rad(cell2mat(bpN_WolfIntactFESNUnC_FAST(:,18:25)));

circstats_median_s = nan(1,8);
circstats_median_m = nan(1,8);
circstats_median_f = nan(1,8);

for i=1:8
    
    circstats_median_s(1,i) = circ_median(legPhaseDiffs_rads_s(:,i));
    circstats_median_m(1,i) = circ_median(legPhaseDiffs_rads_m(:,i));
    circstats_median_f(1,i) = circ_median(legPhaseDiffs_rads_f(:,i));
    
end

% find out max n in bins
for i=1:8
cData = deg2rad(cell2mat(bpN_WolfIntactFESNUnC_SLOW(:,17+i))); % starts at col 18
cData(isnan(cData))=[];
[~,Nslow]=rose(cData,80);
cNslow(i) = max(Nslow);
end

for i=1:8
cData = deg2rad(cell2mat(bpN_WolfIntactFESNUnC_MED(:,17+i))); % starts at col 18
cData(isnan(cData))=[];
[~,Nmed]=rose(cData,80);
cNmed(i) = max(Nmed);
end

for i=1:8
cData = deg2rad(cell2mat(bpN_WolfIntactFESNUnC_FAST(:,17+i))); % starts at col 18
cData(isnan(cData))=[];
[~,Nfast]=rose(cData,80);
cNfast(i) = max(Nfast);
end


% set line colours so goes black, blue...
roseWolfLineCols = [0    0.4470    0.7410
    0    0   0
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330];
set(groot,'defaultAxesColorOrder',roseWolfLineCols);

nameString = {'L1','L2','L3', 'L4', 'R1', 'R2','R3', 'R4'};
plotOrd = [1 3 5 7 2 4 6 8];


f1 = figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad(cell2mat(bpN_WolfIntactFESNUnC_SLOW(:,17+i))); % starts at col 18
cData(isnan(cData))=[];
rose(cData,80);
hold on
c_Median = repelem(circstats_median_s(i),cNslow(i)); % replicates median value to max. no. of bins, so line is same length as maxiumum bin
h=rose(c_Median,500); % change '500' to make wider/narrower
x = get(h,'Xdata');
y = get(h,'Ydata');
g=patch(x,y,'k'); % fills in the median 'bin' in black (k)

title(['Leg phase diff ' nameString{i} '_s']) ;
end

f2 = figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad(cell2mat(bpN_WolfIntactFESNUnC_MED(:,17+i))); % starts at col 18
cData(isnan(cData))=[];
rose(cData,80);
hold on
c_Median = repelem(circstats_median_m(i),cNmed(i)); % replicates median value to max. no. of bins, so line is same length as maxiumum bin
h=rose(c_Median,500); % change '500' to make wider/narrower
x = get(h,'Xdata');
y = get(h,'Ydata');
g=patch(x,y,'k'); % fills in the median 'bin' in black (k)
title(['Leg phase diff ' nameString{i} '_m']) ;
end

f3 = figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad(cell2mat(bpN_WolfIntactFESNUnC_FAST(:,17+i))); % starts at col 18
cData(isnan(cData))=[];
rose(cData,80);
hold on
c_Median = repelem(circstats_median_f(i),cNfast(i)); % replicates median value to max. no. of bins, so line is same length as maxiumum bin
h=rose(c_Median,500); % change '500' to make wider/narrower
x = get(h,'Xdata');
y = get(h,'Ydata');
g=patch(x,y,'k'); % fills in the median 'bin' in black (k)
title(['Leg phase diff ' nameString{i} '_f']) ;
end

