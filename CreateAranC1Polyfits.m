% create fits using polyfit

run PullOutDataForC1Stats.m

%% bpN_data - Aran
[fitresult, gof, output] = createPolyFit(bpN_AranIntact_speed, bpN_AranIntact_dF)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('dutyFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_AranInt_speed-dutyF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_AranInt_speed-dutyF_resid','-dpng','-r0');

close all;

[fitresult, gof, output] = createPolyFit(bpN_AranIntact_speed, bpN_AranIntact_slipF)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('slipFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_AranInt_speed-slipF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_AranInt_speed-slipF_resid','-dpng','-r0');

close all;


[fitresult, gof, output] = createPolyFit(bpN_AranIntact_speed, bpN_AranIntact_strideLen)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('strideLength norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_AranInt_speed-strideLen_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_AranInt_speed-strideLen_resid','-dpng','-r0');

close all;

%% ls_data Aran
%dutyfactor
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_dF)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('dutyFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-dutyF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-dutyF_resid','-dpng','-r0');

close all;

%slipfactor
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_slipF)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('slipFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-slipF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-slipF_resid','-dpng','-r0');

close all;

%stance period
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_stanceP)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('stancePeriod norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-stanceP_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-stanceP_resid','-dpng','-r0');

close all;

%stance X excursion
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_stanceXExc)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('stance X Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-stanceXExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-stanceXExc_resid','-dpng','-r0');

close all;

%stance Y excursion
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_stanceYExc)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('stance Y Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-stanceYExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-stanceYExc_resid','-dpng','-r0');

close all;

%stride period 
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_strideP)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('stridePeriod norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-strideP_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-strideP_resid','-dpng','-r0');

close all;

%swing period 
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_swingP)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('swingPeriod norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-swingP_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-swingP_resid','-dpng','-r0');

close all;

%swing X excursion
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_swingXExc)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('swing X Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-swingXExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-swingXExc_resid','-dpng','-r0');

close all;

%swing Y excursion
[fitresult, gof, output] = createPolyFit(lsN_AranIntact_speed, lsN_AranIntact_swingYExc)
title('fit - Aran Intact');
xlabel('StrideAveVel norm');
ylabel('swing Y Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-swingYExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_AranInt_speed-swingYExc_resid','-dpng','-r0');

close all;
