% calculate correlation with speed for circular data - wolf uncooled
run PullOutDataForC1Stats.m

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

%bp_data
[r1, p1] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_StrideAveVelAng), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_StrideAveVelAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride mean velocity angle (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-StrideAveVelAng_scat','-dpng','-r0');
close all

[r2, p2] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_StrideDeltaVelAng), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_StrideDeltaVelAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride delta velocity angle (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-StrideDeltaVelAng_scat','-dpng','-r0');
close all

[r3, p3] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_StrideAveYawAng), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_StrideAveYawAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride mean yaw angle (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-StrideAveYawAng_scat','-dpng','-r0');
close all

[r4, p4] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_StrideDeltaYawAng), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_StrideDeltaYawAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride delta yaw angle (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-StrideDeltaYawAng_scat','-dpng','-r0');
close all

[r5, p5] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffL1), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffL1,20,SpidColors(1,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference L1 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffL1_scat','-dpng','-r0');
close all

[r6, p6] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffL2), bpN_WolfIntactFESNUnC_vel);

[r7, p7] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffL3), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffL3,20,SpidColors(3,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference L3 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffL3_scat','-dpng','-r0');
close all

[r8, p8] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffL4), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffL4,20,SpidColors(4,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference L4 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffL4_scat','-dpng','-r0');
close all

[r9, p9] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffR1), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffR1,20,SpidColors(5,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference R1 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffR1_scat','-dpng','-r0');
close all

[r10, p10] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffR2), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffR2,20,SpidColors(6,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference R2 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffR2_scat','-dpng','-r0');
close all

[r11, p11] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffR3), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffR3,20,SpidColors(7,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference R3 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffR3_scat','-dpng','-r0');
close all

[r12, p12] = circ_corrcl(deg2rad(c_bpN_WolfIntactFESNUnC_legPhaseDiffR4), bpN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(bpN_WolfIntactFESNUnC_vel,c_bpN_WolfIntactFESNUnC_legPhaseDiffR4,20,SpidColors(8,:),'filled');
xlabel('Stride mean velocity (norm)');
ylabel('Leg phase difference R4 (degrees)');
title('Wolf intact uncooled (bpN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('bp_WolfIntFESNUnC_vel-legPhaseDiffR4_scat','-dpng','-r0');
close all

% create column of variable names
rP_bp_varColumn = {'StrideAveVelAng';'StrideDeltaVelAng';'StrideAveYawAng';'StrideDeltaYawAng';'legPhaseDiffL1';'legPhaseDiffL2';'legPhaseDiffL3';...
    'legPhaseDiffL4';'legPhaseDiffR1';'legPhaseDiffR2';'legPhaseDiffR3';'legPhaseDiffR4'};

% create header row
rP_bp_dataHeaders = {'variables','rho','pval'};

% put values together - each row = 1 variable
rP_bp_values = {r1, p1; r2, p2; r3, p3; r4, p4; r5, p5; r6, p6; r7, p7; r8, p8; r9, p9; r10, p10; r11, p11; r12, p12};

% % Here I tried to edit the for loop to add a column of significance but
% it didn't really work (and I don't really need it)
% rP_bp_values_temp = cell2mat(rP_bp_values);
% sigs = [];
% 
% for i=1:length(rP_bp_values_temp);
%     if rP_bp_values_temp(i,2) < 0.05
%         sigs(i) = 1;
%     else
%         sigs(i) = 0;
%     end
% end
% 
% sigs = num2cell(sigs);
% sigs = sigs';

% rP_bpN_WolfIntactFESNUnC = cell(1,length(rsPs_bpN_WolfIntactFESNUnC));
% [rP_bpN_WolfIntactFESNUnC{1,1:end}] = deal(rP_bp_headers{:});
% rP_bpN_WolfIntactFESNUnC = vertcat(rP_bpN_WolfIntactFESNUnC,num2cell(rsPs_bpN_WolfIntactFESNUnC),num2cell(sigs));


% put together left had variable column with data values
rsPs_bpN_WolfIntactFESNUnC_temp = horzcat(rP_bp_varColumn,rP_bp_values);

% add header row
rsPs_bpN_WolfIntactFESNUnC = vertcat(rP_bp_dataHeaders,rsPs_bpN_WolfIntactFESNUnC_temp);


cell2csv([saveDir,'speedEffectData_bpN_WolfIntactFESNUnC.csv'], rsPs_bpN_WolfIntactFESNUnC);

clear r* p* rsPs_bpN_WolfIntactFESNUnC_temp rP_bp_values_temp

%ls_data
[r1, p1] = circ_corrcl(deg2rad(c_lsN_WolfIntactFESNUnC_StrideAveVelAng), lsN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(lsN_WolfIntactFESNUnC_vel,c_lsN_WolfIntactFESNUnC_StrideAveVelAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride mean velocity angle (degrees)');
title('Wolf intact uncooled (lsN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfIntFESNUnC_vel-StrideAveVelAng_scat','-dpng','-r0');
close all

[r2, p2] = circ_corrcl(deg2rad(c_lsN_WolfIntactFESNUnC_StrideDeltaVelAng), lsN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(lsN_WolfIntactFESNUnC_vel,c_lsN_WolfIntactFESNUnC_StrideDeltaVelAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride delta velocity angle (degrees)');
title('Wolf intact uncooled (lsN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfIntFESNUnC_vel-StrideDeltaVelAng_scat','-dpng','-r0');
close all

[r3, p3] = circ_corrcl(deg2rad(c_lsN_WolfIntactFESNUnC_StrideAveYawAng), lsN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(lsN_WolfIntactFESNUnC_vel,c_lsN_WolfIntactFESNUnC_StrideAveYawAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride mean yaw angle (degrees)');
title('Wolf intact uncooled (lsN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfIntFESNUnC_vel-StrideAveYawAng_scat','-dpng','-r0');
close all

[r4, p4] = circ_corrcl(deg2rad(c_lsN_WolfIntactFESNUnC_StrideDeltaYawAng), lsN_WolfIntactFESNUnC_vel);
f1 = figure();
set(f1,'units', 'normalized'); set(f1,'Position', [0 0.0364583 1 0.875]);
f1 = scatter(lsN_WolfIntactFESNUnC_vel,c_lsN_WolfIntactFESNUnC_StrideDeltaYawAng,10,'b','filled');
xlabel('Stride mean velocity (norm)');
ylabel('Stride delta yaw angle (degrees)');
title('Wolf intact uncooled (lsN)');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('ls_WolfIntFESNUnC_vel-StrideDeltaYawAng_scat','-dpng','-r0');
close all


rP_ls_varColumn = {'StrideAveVelAng';'StrideDeltaVelAng';'StrideAveYawAng';'StrideDeltaYawAng'};

% create header row
rP_ls_dataHeaders = {'variables','rho','pval'};

% put values together - each row = 1 variable
rp_ls_values = {r1, p1; r2, p2; r3, p3; r4, p4};

% put together left had variable column with data values
rsPs_lsN_WolfIntactFESNUnC_temp = horzcat(rP_ls_varColumn,rp_ls_values);

% add header row
rsPs_lsN_WolfIntactFESNUnC = vertcat(rP_ls_dataHeaders,rsPs_lsN_WolfIntactFESNUnC_temp);


cell2csv([saveDir,'speedEffectData_lsN_WolfIntactFESNUnC.csv'], rsPs_lsN_WolfIntactFESNUnC);

clear r* p* rP_ls_varColumn rsPs_lsN_WolfIntactFESNUnC_temp rP_ls_dataHeaders


