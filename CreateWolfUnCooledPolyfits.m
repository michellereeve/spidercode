% create fits using polyfit

run PullOutDataForC1Stats.m

saveDir = ['/Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/TrendsWithSpeed/NonCirc_TrendsWithSpeed/WolfPolyfits_nonCirc_UnC_T3/StatsTables/'];
cd /Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/TrendsWithSpeed/NonCirc_TrendsWithSpeed/WolfPolyfits_nonCirc_UnC_T3/
% Just export ls data

%% ls_data Wolf uncooled

% create empty cell array for values to go into
r2_rho_p_ls_values = cell(9,3);
% f-statistic p-value is from regstats -> stats.fstat.pval
% slope/gradient is from fitresult -> fitresult.p1 (can check this by doing formula(fitresult) -> gives p1*x + p2
% r-squared value from gof.rsquare

%dutyfactor
stats = regstats(lsN_WolfIntactFESNUnC_dF,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{1,3} = fstatpval; % p-value

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_dF,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('dutyFactor');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-dutyF_fit','-dpng','-r0');
r2_rho_p_ls_values{1,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{1,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-dutyF_resid','-dpng','-r0');

clear stats;
close all;


%slipfactor
stats = regstats(lsN_WolfIntactFESNUnC_slipF,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{2,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_slipF,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('slipFactor');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-slipF_fit','-dpng','-r0');

r2_rho_p_ls_values{2,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{2,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-slipF_resid','-dpng','-r0');

clear stats;
close all;


%stance period
stats = regstats(lsN_WolfIntactFESNUnC_stanceP,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{3,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceP,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('stancePeriod (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-stanceP_fit','-dpng','-r0');

r2_rho_p_ls_values{3,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{3,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-stanceP_resid','-dpng','-r0');

clear stats;
close all;


%stance X excursion
stats = regstats(lsN_WolfIntactFESNUnC_stanceXExc,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{4,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceXExc,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('stance X Excursion (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-stanceXExc_fit','-dpng','-r0');

r2_rho_p_ls_values{4,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{4,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-stanceXExc_resid','-dpng','-r0');

clear stats;
close all;


%stance Y excursion
stats = regstats(lsN_WolfIntactFESNUnC_stanceYExc,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{5,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_stanceYExc,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('stance Y Excursion (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-stanceYExc_fit','-dpng','-r0');

r2_rho_p_ls_values{5,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{5,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-stanceYExc_resid','-dpng','-r0');

clear stats;
close all;


%stride period 
stats = regstats(lsN_WolfIntactFESNUnC_strideP,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{6,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_strideP,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('stridePeriod (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-strideP_fit','-dpng','-r0');

r2_rho_p_ls_values{6,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{6,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-strideP_resid','-dpng','-r0');

clear stats;
close all;


%swing period 
stats = regstats(lsN_WolfIntactFESNUnC_swingP,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{7,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingP,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('swingPeriod (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-swingP_fit','-dpng','-r0');

r2_rho_p_ls_values{7,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{7,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-swingP_resid','-dpng','-r0');

clear stats;
close all;


%swing X excursion
stats = regstats(lsN_WolfIntactFESNUnC_swingXExc,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{8,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingXExc,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('swing X Excursion (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-swingXExc_fit','-dpng','-r0');

r2_rho_p_ls_values{8,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{8,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-swingXExc_resid','-dpng','-r0');

clear stats;
close all;


%swing Y excursion
stats = regstats(lsN_WolfIntactFESNUnC_swingYExc,lsN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

r2_rho_p_ls_values{9,3} = fstatpval;

[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFESNUnC_vel, lsN_WolfIntactFESNUnC_swingYExc,fstatpval);
% title('fit - Wolf uncooled Intact');
% xlabel('stride mean velocity (norm)');
% ylabel('swing Y Excursion (norm)');
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-swingYExc_fit','-dpng','-r0');

r2_rho_p_ls_values{9,2} = fitresult.p1; % gradient
r2_rho_p_ls_values{9,1} = gof.rsquare;

% % plot residuals with line at y=0
% f2 = figure('Name', 'residuals');
% set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
% plot(output.residuals,'-ko')
% legend('residuals');
% refline(0)
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('ls_WolfIntUnC_vel-swingYExc_resid','-dpng','-r0');

clear stats;
close all;

% create column of variable names
r2_rho_p_ls_varColumn = {'duty factor';'slip factor';'stance period';'stance fore-aft exc.';'stance med-lat exc.';'stride period';'swing period';'swing fore-aft exc.';...
    'swing med-lat exc.'};
% create header row
r2_rho_p_ls_dataHeaders = {'variables','R^2','slope','p-value'};

% put together left had variable column with data values
r2_rho_p_lsN_WolfIntactFESNUnC_temp = horzcat(r2_rho_p_ls_varColumn,r2_rho_p_ls_values);

% add header row
r2_rho_p_lsN_WolfIntactFESNUnC = vertcat(r2_rho_p_ls_dataHeaders,r2_rho_p_lsN_WolfIntactFESNUnC_temp);

% save to csv
cell2csv([saveDir,'speedEffectData_noncirc_lsN_WolfIntactFESNUnC.csv'], r2_rho_p_lsN_WolfIntactFESNUnC);


%% bpN_data - Wolf Uncooled

% dutyfactor
stats = regstats(bpN_WolfIntactFESNUnC_dF,bpN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult, gof, output] = createPolyFit(bpN_WolfIntactFESNUnC_vel, bpN_WolfIntactFESNUnC_dF,fstatpval)
title('fit - Wolf uncooled Intact');
xlabel('stride mean velocity (norm)');
ylabel('dutyFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntUnC_vel-dutyF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntUnC_vel-dutyF_resid','-dpng','-r0');

clear stats;
close all;

%slipFactor
stats = regstats(bpN_WolfIntactFESNUnC_slipF,bpN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult, gof, output] = createPolyFit(bpN_WolfIntactFESNUnC_vel, bpN_WolfIntactFESNUnC_slipF,fstatpval)
title('fit - Wolf uncooled Intact');
xlabel('stride mean velocity (norm)');
ylabel('slipFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntUnC_vel-slipF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntUnC_vel-slipF_resid','-dpng','-r0');

clear stats;
close all;

%stride length
stats = regstats(bpN_WolfIntactFESNUnC_strideLen,bpN_WolfIntactFESNUnC_vel,'linear');
fstatpval = stats.fstat.pval;

[fitresult, gof, output] = createPolyFit(bpN_WolfIntactFESNUnC_vel, bpN_WolfIntactFESNUnC_strideLen,fstatpval)
title('fit - Wolf uncooled Intact');
xlabel('stride mean velocity (norm)');
ylabel('strideLength (norm)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntUnC_vel-strideLen_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntUnC_vel-strideLen_resid','-dpng','-r0');

clear stats;
close all;

