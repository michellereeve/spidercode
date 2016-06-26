% ANOVA to test the effect of segNum, lVsR, legNum and individual, wolf flat intact uncooled - MED speed bin
% USED JUNE 2016

run PullOutDataForC1Stats.m
clear lsN_is* bpN_is* lsN_Aran* bpN_Aran* c_bpN_Aran* c_lsN_Aran* *FESNC_* *FESNC % clear all non- uncooled stuff and all aran
% group into speed bins
[lsN_WolfIntactFESNUnC_SLOW, lsN_WolfIntactFESNUnC_MED, lsN_WolfIntactFESNUnC_FAST] = SetSpeedBins(lsN_WolfIntactFESNUnC,lsN_WolfIntactFESNUnC_vel);

% FOR SLOW SPEED BIN
stanceXExc_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,19)); % col 19 = stanceXExcur
stanceYExc_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,21));
swingXExc_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,20)); 
swingYExc_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,22));
dutyFactor_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,18));
slipFactor_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,23));
swingPeriod_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,16));
stancePeriod_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,17));

%legNumEffect_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,7)); % col 7 = legNum
indivEffect_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,3)); % col 3 - subjectNum
segNumEffect_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,24)); % col 24 = segNum
%LvREffect_S = cell2mat(lsN_WolfIntactFESNUnC_SLOW(:,25)); % col 25 = RvL

% FOR MED SPEED BIN
stanceXExc_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,19)); % col 19 = stanceXExcur
stanceYExc_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,21));
swingXExc_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,20)); 
swingYExc_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,22));
dutyFactor_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,18));
slipFactor_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,23));
swingPeriod_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,16));
stancePeriod_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,17));

%legNumEffect_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,7)); % col 7 = legNum
indivEffect_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,3)); % col 3 - subjectNum
segNumEffect_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,24)); % col 24 = segNum
%LvREffect_M = cell2mat(lsN_WolfIntactFESNUnC_MED(:,25)); % col 25 = RvL

% FOR FAST SPEED BIN
stanceXExc_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,19)); % col 19 = stanceXExcur
stanceYExc_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,21));
swingXExc_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,20)); 
swingYExc_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,22));
dutyFactor_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,18));
slipFactor_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,23));
swingPeriod_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,16));
stancePeriod_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,17));

%legNumEffect_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,7)); % col 7 = legNum
indivEffect_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,3)); % col 3 - subjectNum
segNumEffect_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,24)); % col 24 = segNum
%LvREffect_F = cell2mat(lsN_WolfIntactFESNUnC_FAST(:,25)); % col 25 = RvL


% % ANOVA with just LegNum
% [pVals, table, stats] = anovan(stanceXExc,{legNumEffect}, 'model',2, 'varnames', {'legNum'});
% %[pVals, table, stats] = anova1(stanceXExc,legNumEffect);
% 
% %ANOVA with LegNum & individual
% [pVals, table, stats] = anovan(stanceXExc,{legNumEffect,indivEffect}, 'model',2, 'random', 2, 'varnames', {'legNum' 'indiv'});
% % NOT FULL RANK error - so trying with the below instead

% % ANOVA with segNum, LvsR, and individual
% [pVals, table, stats] = anovan(stanceXExc,{segNumEffect,LvREffect, indivEffect}, 'model',2, 'random', 3, 'varnames', {'segNum' 'LvsR' 'indiv'});
% 
% % as above for only interaction of interest
% effects = [segNumEffect,LvREffect,indivEffect];
% interactionCoding = [1 1 0; 1 0 0; 0 1 0; 0 0 1];
% 
% % for STANCE X
% [pVals, table, stats] = anovan(stanceXExc,effects,'random', [3], 'model', interactionCoding,'varnames', {'segNum' 'LvsR' 'indiv'});
% 
%  %post hoc analyses
%     figure;
%     [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
%     figure(hCont);     
%     title(['Click on group to test: stanceXExc']);
%   
%     
% % FOR STANCE Y 
%     [pVals, table, stats] = anovan(stanceYExc,effects,'random', [3], 'model', interactionCoding,'varnames', {'segNum' 'LvsR' 'indiv'});
% 
%  %post hoc analyses
%     figure;
%     [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
%     figure(hCont);     
%     title(['Click on group to test: stanceYExc']);
%     
% % % aNOVA for segNum + ind (no interacton)
% % [pVals, table, stats] = anovan(stanceXExc,{segNumEffect,indivEffect}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'});
% % 
% %  %post hoc analyses
% %     figure;
% %     [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1 2]); %#ok<ASGLU>
% %     figure(hCont);     
% %     title(['Click on group to test: stanceXExc']);
%     
% % FOR SWING X
%     [pVals, table, stats] = anovan(swingXExc,effects,'random', [3], 'model', interactionCoding,'varnames', {'segNum' 'LvsR' 'indiv'});
% 
%  %post hoc analyses
%     figure;
%     [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1 2]); %#ok<ASGLU>
%     figure(hCont);     
%     title(['Click on group to test: swingXExc']);
%     
% % FOR SWING Y
%     [pVals, table, stats] = anovan(swingYExc,effects,'random', [3], 'model', interactionCoding,'varnames', {'segNum' 'LvsR' 'indiv'});
% 
%  %post hoc analyses
%     figure;
%     [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
%     figure(hCont);     
%     title(['Click on group to test: swingYExc']);    

% FINAL STUFF - JUST SEGNUM & INDIVIDUAL FOR WOLF INTACT NON-CIRC DATA - NO INTERACTION

%% MED DATA
% STANCE X
[pVals, table, stats] = anovan(stanceXExc_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceXExc']);    

% STANCE Y    
[pVals, table, stats] = anovan(stanceYExc_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceYExc']);  

% SWING X    
[pVals, table, stats] = anovan(swingXExc_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: swingXExc']);

% SWING Y    
[pVals, table, stats] = anovan(swingYExc_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: swingYExc']);

% DUTY FACTOR    
[pVals, table, stats] = anovan(dutyFactor_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: dutyFactor']);

% SLIP FACTOR    
[pVals, table, stats] = anovan(slipFactor_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 
% no sig results

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: slipFactor']);

%SWING PERIOD   
[pVals, table, stats] = anovan(swingPeriod_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% fH = boxplot(swingPeriod, segNumEffect); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: swingPeriod']);

% STANCE PERIOD
[pVals, table, stats] = anovan(stancePeriod_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 


figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stancePeriod']);
    

% useful code for converting data for anova 
%[groupIndexArray] = grp2idx(groupCellArray(:,1));


%% SLOW DATA
% STANCE X
[pVals, table, stats] = anovan(stanceXExc_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: stanceXExc']);    

% STANCE Y    
[pVals, table, stats] = anovan(stanceYExc_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: stanceYExc']);  

% SWING X    
[pVals, table, stats] = anovan(swingXExc_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: swingXExc']);

% SWING Y    
[pVals, table, stats] = anovan(swingYExc_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: swingYExc']);

% DUTY FACTOR    
[pVals, table, stats] = anovan(dutyFactor_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: dutyFactor']);

% SLIP FACTOR    
[pVals, table, stats] = anovan(slipFactor_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: slipFactor']);

%SWING PERIOD   
[pVals, table, stats] = anovan(swingPeriod_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: swingPeriod']);

% STANCE PERIOD
[pVals, table, stats] = anovan(stancePeriod_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: stancePeriod'])

%% FAST DATA
% STANCE X
[pVals, table, stats] = anovan(stanceXExc_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceXExc']);    

% STANCE Y    
[pVals, table, stats] = anovan(stanceYExc_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceYExc']);  

% SWING X    
[pVals, table, stats] = anovan(swingXExc_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: swingXExc']);

% SWING Y    
[pVals, table, stats] = anovan(swingYExc_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: swingYExc']);

% DUTY FACTOR    
[pVals, table, stats] = anovan(dutyFactor_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: dutyFactor']);

% SLIP FACTOR    
[pVals, table, stats] = anovan(slipFactor_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 
% no sig results

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: slipFactor']);

%SWING PERIOD   
[pVals, table, stats] = anovan(swingPeriod_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

% figure;
% [cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
% figure(hCont);     
% title(['Click on group to test: swingPeriod']);

% STANCE PERIOD
[pVals, table, stats] = anovan(stancePeriod_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stancePeriod'])

%% PRETTY MULTCOMPARE PLOTS

% SLOW SLIP FACTOR
[pVals, table, stats] = anovan(slipFactor_S,{segNumEffect_S,indivEffect_S}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 
figure;
[cCont,mCont,hCont,gnamesCont] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: slipFactor']);

dat_s = flipud(mCont(:,1)); % flipud reverses order - to make plots be front = right
err_s = flipud(mCont(:,2));
f1=figure;
errorbar(dat_s,err_s,'bx');
ax=gca;
ax.XTick = [1 2 3 4]
box(ax,'off');
set(gca,'xticklabel',flipud(gnamesCont))
%title('SLOW slip factor');
ylabel('slip factor');


% STANCE X MED vs FAST
[pVals, table, stats] = anovan(stanceXExc_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 
figure;
[cCont,mCont_m,hCont,gnamesCont_m] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceXExc']); 

[pVals, table, stats] = anovan(stanceXExc_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 
figure;
[cCont,mCont_f,hCont,gnamesCont_f] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceXExc']); 

dat_m = flipud(mCont_m(:,1));
err_m = flipud(mCont_m(:,2));
dat_f = flipud(mCont_f(:,1));
err_f = flipud(mCont_f(:,2));

figure;
sp1=subplot(1,2,1);
errorbar(dat_m,err_m,'bx'); % be nice to do points 1 col, errbars another & weights
ylabel('stance fore-aft excursion (norm.)');
ax=gca;
box(ax,'off');
ax.XTick = [1 2 3 4]
set(gca,'xticklabel',flipud(gnamesCont_m));
title('MEDIUM');

sp2=subplot(1,2,2);
xlim([0.85 1.40]);
errorbar(dat_f,err_f,'bx');
ax=gca;
box(ax,'off');
ax.XTick = [1 2 3 4]
ax.YTickLabel = '';
set(gca,'xticklabel',flipud(gnamesCont_f));
title('FAST');
linkaxes([sp2,sp1],'y'); % links subplot axes so they show same scale on y axis, based on sp2 data (this is first)


% STANCE Y    
[pVals, table, stats] = anovan(stanceYExc_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont_m,hCont,gnamesCont_m] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceYExc']); 

[pVals, table, stats] = anovan(stanceYExc_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont_f,hCont,gnamesCont_f] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: stanceYExc']);  

dat_m = flipud(mCont_m(:,1));
err_m = flipud(mCont_m(:,2));
dat_f = flipud(mCont_f(:,1));
err_f = flipud(mCont_f(:,2));
figure;
sp1=subplot(1,2,1);
errorbar(dat_m,err_m,'bx'); % be nice to do points 1 col, errbars another & weights
ylabel('stance med-lat excursion (norm.)');  % doesn't work - why?
ax=gca;
box(ax,'off');
ax.XTick = [1 2 3 4];
set(gca,'xticklabel',flipud(gnamesCont_m));
title('MEDIUM');
sp2=subplot(1,2,2);
errorbar(dat_f,err_f,'bx');
ax=gca;
box(ax,'off');
ax.XTick = [1 2 3 4]
ax.YTickLabel = '';
set(gca,'xticklabel',flipud(gnamesCont_f));
title('FAST');
linkaxes([sp2,sp1],'y');

% DutyFactor
[pVals, table, stats] = anovan(dutyFactor_M,{segNumEffect_M,indivEffect_M}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont_m,hCont,gnamesCont_m] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: dutyFactor']);
 
[pVals, table, stats] = anovan(dutyFactor_F,{segNumEffect_F,indivEffect_F}, 'model',1, 'random', 2, 'varnames', {'segNum' 'indiv'}); 

figure;
[cCont,mCont_f,hCont,gnamesCont_f] = multcompare(stats,'ctype', 'bonferroni','display','on','dimension',[1]); %#ok<ASGLU>
figure(hCont);     
title(['Click on group to test: dutyFactor']);

dat_m = flipud(mCont_m(:,1));
err_m = flipud(mCont_m(:,2));
dat_f = flipud(mCont_f(:,1));
err_f = flipud(mCont_f(:,2));
figure;
sp1=subplot(1,2,1);
errorbar(dat_m,err_m,'bx'); % be nice to do points 1 col, errbars another & weights
ylabel('duty factor');  % doesn't work - why?
ax=gca;
box(ax,'off');
ax.XTick = [1 2 3 4]
ylim([0.24 0.40]);
set(gca,'xticklabel',flipud(gnamesCont_m));
title('MEDIUM');
sp2=subplot(1,2,2);
errorbar(dat_f,err_f,'bx');
ax=gca;
box(ax,'off');
ax.XTick = [1 2 3 4]
ylim([0.24 0.40]);
ax.YTickLabel = '';
set(gca,'xticklabel',gnamesCont_f)
title('FAST');
%linkaxes([sp1,sp2],'y'); didnt work for this one
