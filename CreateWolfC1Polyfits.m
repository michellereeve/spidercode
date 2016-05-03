% create fits using polyfit

run PullOutDataForC1Stats.m

%% bpN_data - Wolf
[fitresult, gof, output] = createPolyFit(bpN_WolfIntactFlatESN_speed, bpN_WolfIntactFlatESN_dF)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('dutyFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfInt_speed-dutyF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfInt_speed-dutyF_resid','-dpng','-r0');

close all;

[fitresult, gof, output] = createPolyFit(bpN_WolfIntactFlatESN_speed, bpN_WolfIntactFlatESN_slipF)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('slipFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfInt_speed-slipF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfInt_speed-slipF_resid','-dpng','-r0');

close all;


[fitresult, gof, output] = createPolyFit(bpN_WolfIntactFlatESN_speed, bpN_WolfIntactFlatESN_strideLen)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('strideLength norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfInt_speed-strideLen_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfInt_speed-strideLen_resid','-dpng','-r0');

close all;

%% ls_data Wolf
%dutyfactor
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_dF)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('dutyFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-dutyF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-dutyF_resid','-dpng','-r0');

close all;

%slipfactor
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_slipF)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('slipFactor');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-slipF_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-slipF_resid','-dpng','-r0');

close all;

%stance period
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_stanceP)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('stancePeriod norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-stanceP_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-stanceP_resid','-dpng','-r0');

close all;

%stance X excursion
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_stanceXExc)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('stance X Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-stanceXExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-stanceXExc_resid','-dpng','-r0');

close all;

%stance Y excursion
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_stanceYExc)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('stance Y Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-stanceYExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-stanceYExc_resid','-dpng','-r0');

close all;

%stride period 
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_strideP)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('stridePeriod norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-strideP_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-strideP_resid','-dpng','-r0');

close all;

%swing period 
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_swingP)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('swingPeriod norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-swingP_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-swingP_resid','-dpng','-r0');

close all;

%swing X excursion
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_swingXExc)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('swing X Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-swingXExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfInt_speed-swingXExc_resid','-dpng','-r0');

close all;

%swing Y excursion
[fitresult, gof, output] = createPolyFit(lsN_WolfIntactFlatESN_speed, lsN_WolfIntactFlatESN_swingYExc)
title('fit - Wolf Intact');
xlabel('StrideAveVel norm');
ylabel('swing Y Excursion norm');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfIntFESN_speed-swingYExc_fit','-dpng','-r0');

% plot residuals with line at y=0
f2 = figure('Name', 'residuals');
set(f2,'units', 'normalized'); set(f2,'Position', [0 0.0364583 1 0.875]);
plot(output.residuals,'-ko')
legend('residuals');
refline(0)
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfIntFESN_speed-swingYExc_resid','-dpng','-r0');

close all;
