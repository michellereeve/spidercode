% create multipanel plot for sig results

% originally had 8 plots: stride length, dF, stanceP, stance fore-aft,
% strideP, swingP, swing fore-aft, swing med-lat 
% changing to just interesting sig results 15062016: stanceP, swingP, dF,
% stance fore-aft

run PullOutDataForC1Stats.m
%[lsN_WolfIntactFESNUnC_SLOW, lsN_WolfIntactFESNUnC_MED, lsN_WolfIntactFESNUnC_FAST] = SetSpeedBins(lsN_WolfIntactFESNUnC,lsN_WolfIntactFESNUnC_vel);

%% NON CIRCULAR DATA
cd /Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/TrendsWithSpeed/NonCirc_TrendsWithSpeed/WolfPolyfits_nonCirc_UnC_T3/
f1 = figure;

%stance period
stats = regstats(lsN_WolfIntactFESNUnC_stanceP,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult3, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceP,fstatpval)

figure(f1);
sp1=subplot(4,1,1);
plot(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceP,'.k');
hold on
plot(fitresult3,'-b');
xlabel('');
ylabel('stance period (norm)');
legend('data','fitted curve');
legend(sp1,'boxoff');
xlim([0 2.52]);
box(sp1,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';


%swing period
stats = regstats(lsN_WolfIntactFESNUnC_swingP,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult6, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingP,fstatpval)

figure(f1);
sp2=subplot(4,1,2);
plot(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingP,'.k'); % add 'predobs' argument if want prediction bounds
hold on
plot(fitresult6,'r');
xlabel('');
ylabel('swing period (norm)');
legend(sp2,'hide');
xlim([0 2.52]);
box(sp2,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%dutyfactor
stats = regstats(lsN_WolfIntactFESNUnC_dF,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult2, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_dF,fstatpval)
figure(f1);
sp3=subplot(4,1,3);
plot(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_dF,'.b');
hold on
plot(fitresult2,'k');
xlabel('');
ylabel('duty factor');
legend(sp3,'hide');
xlim([0 2.52]);
box(sp3,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%stance X excur
stats = regstats(lsN_WolfIntactFESNUnC_stanceXExc,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult4, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceXExc,fstatpval)

figure(f1);
sp4=subplot(4,1,4);
plot(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceXExc,'.k');
hold on
plot(fitresult4,'b');
xlabel('Stride mean velocity (norm)');
ylabel('stance fore-aft excur. (norm)');
legend(sp4,'hide');
xlim([0 2.52]);
box(sp4,'off');
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.005 0.035];


% %swing X excur
% stats = regstats(lsN_WolfIntactFESNUnC_swingXExc,lsN_WolfIntactFESNUnC_vel,'linear');
% fstatpval = stats.fstat.pval;
% 
% [fitresult7, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingXExc,fstatpval)
% 
% figure(f1);
% subplot(4,2,7);
% plot(fitresult7, lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingXExc,'predobs');
% xlabel('stride mean velocity (norm)');
% ylabel('swing fore-aft excur. (norm)');
% 
% %swing Y excur
% stats = regstats(lsN_WolfIntactFESNUnC_swingXExc,lsN_WolfIntactFESNUnC_vel,'linear');
% fstatpval = stats.fstat.pval;
% 
% [fitresult8, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingYExc,fstatpval)
% 
% figure(f1);
% subplot(4,2,8);
% plot(fitresult8, lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingYExc,'predobs');
% xlabel('stride mean velocity (norm)');
% ylabel('swing med-lat excur. (norm)');

% %stride length
% stats = regstats(bpN_WolfIntactFESNUnC_strideLen,bpN_WolfIntactFESNUnC_vel,'linear');
% fstatpval = stats.fstat.pval;
% 
% [fitresult1, gof, output] = createPolyFit(bpN_WolfIntactFESNUnC_vel, bpN_WolfIntactFESNUnC_strideLen,fstatpval)
% 
% f1 = figure;
% figure(f1);
% subplot(4,2,1);
% plot(fitresult1, bpN_WolfIntactFESNUnC_vel, bpN_WolfIntactFESNUnC_strideLen,'predobs');
% xlabel('');
% ylabel('strideLength (norm)');

% %stride period
% stats = regstats(lsN_WolfIntactFESNUnC_strideP,lsN_WolfIntactFESNUnC_vel,'linear');
% fstatpval = stats.fstat.pval;
% 
% [fitresult5, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_strideP,fstatpval)
% 
% figure(f1);
% subplot(4,2,5);
% plot(fitresult5, lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_strideP,'predobs');
% xlabel('');
% ylabel('stride period (norm)');

%% CIRCULAR DATA

saveDir = ['/Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/TrendsWithSpeed/Circ_TrendsWithSpeed/StatsTables/'];
cd /Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/TrendsWithSpeed/Circ_TrendsWithSpeed/

SpidColors = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0 0 0.4828];

f2 = figure;

% stride ave vel ang
figure(f2);
sp01=subplot(5,2,1);
plot(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_StrideAveVelAng,'.k');
xlabel('');
ylabel('stride mean vel. ang. (deg)');
%xlim([0 2.52]);
box(sp01,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';


% stride delta vel ang
figure(f2);
sp02=subplot(5,2,2);
plot(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_StrideDeltaVelAng,'.k');
xlabel('');
ylabel('stride delta vel. ang. (deg)');
%xlim([0 2.52]);
box(sp02,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';


%leg phase diff L1
figure(f2);
sp03=subplot(5,2,3);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffL1,10,SpidColors(1,:),'filled');
xlabel('');
ylabel('leg phase diff. L1 (deg)');
%xlim([0 2.52]);
box(sp03,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%leg phase diff L3
figure(f2);
sp04=subplot(5,2,7);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffL3,10,SpidColors(3,:),'filled');
xlabel('');
ylabel('leg phase diff. L3 (deg)');
%xlim([0 2.52]);
box(sp04,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%leg phase diff L4
figure(f2);
sp05=subplot(5,2,9);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffL3,10,SpidColors(4,:),'filled');
xlabel('stride mean velocity (norm)');
ylabel('leg phase diff. L4 (deg)');
%xlim([0 2.52]);
box(sp05,'off');
ax = gca;
%ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%leg phase diff R1
figure(f2);
sp06=subplot(5,2,4);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffR1,10,SpidColors(5,:),'filled');
xlabel('');
ylabel('leg phase diff. R1 (deg)');
%xlim([0 2.52]);
box(sp06,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%leg phase diff R2
figure(f2);
sp07=subplot(5,2,6);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffR2,10,SpidColors(6,:),'filled');
xlabel('');
ylabel('leg phase diff. R2 (deg)');
%xlim([0 2.52]);
box(sp07,'off');
ax = gca;
ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%leg phase diff R3
figure(f2);
sp08=subplot(5,2,8);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffR3,10,SpidColors(7,:),'filled');
xlabel('');
ylabel('leg phase diff. R3 (deg)');
%xlim([0 2.52]);
box(sp08,'off');
ax = gca;
%ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';

%leg phase diff R4
figure(f2);
sp09=subplot(5,2,10);
scatter(bpN_WolfIntactFESNUnC_vel, c_bpN_WolfIntactFESNUnC_legPhaseDiffR4,10,SpidColors(8,:),'filled');
xlabel('stride mean velocity (norm)');
ylabel('leg phase diff. R4 (deg)');
%xlim([0 2.52]);
box(sp09,'off');
ax = gca;
%ax.XTickLabel = '';
ax.TickLength = [0.005 0.035];
ax.TickDir = 'out';