%% Wolf intact flat ESN uncooled - SLOW
run PullOutDataForC1Stats.m

saveDir = ['/Users/Michelle/Research/data/__DATA_ANALYSIS/Statistics/'];
clear lsN_is* bpN_is* lsN_Aran* bpN_Aran* c_bpN_Aran* c_lsN_Aran* *FESNC_* *FESNC % clear all non- uncooled stuff and all aran 

[bpN_WolfIntactFESNUnC_SLOW, bpN_WolfIntactFESNUnC_MED, bpN_WolfIntactFESNUnC_FAST] = SetSpeedBins(bpN_WolfIntactFESNUnC,bpN_WolfIntactFESNUnC_vel);

%% CALCULATE MEDIANS + SD
% pre-allocate variables
circstats_all = nan(24,3);
circstats_median = nan(24,1);
circstats_std = nan(24,1);
circstats_speeds = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3]; % 1=slow 2=med 3=fast

% also create separate arrays
circstats_slow = nan(8,2);
circstats_med = nan(8,2);
circstats_fast = nan(8,2);

legPhaseDiffs_rads_s = circ_ang2rad(cell2mat(bpN_WolfIntactFESNUnC_SLOW(:,18:25)));
legPhaseDiffs_rads_m = circ_ang2rad(cell2mat(bpN_WolfIntactFESNUnC_MED(:,18:25)));
legPhaseDiffs_rads_f = circ_ang2rad(cell2mat(bpN_WolfIntactFESNUnC_FAST(:,18:25)));

for i=1:8
    
    circstats_median(i,1) = mod(circ_rad2ang(circ_median(legPhaseDiffs_rads_s(:,i))),360); % mod 360 to convert negative angles to eqiavalent between 0-360
    circstats_std(i,1) = circ_rad2ang(circ_std(legPhaseDiffs_rads_s(:,i)));
    
    circstats_slow(i,1) = mod(circ_rad2ang(circ_median(legPhaseDiffs_rads_s(:,i))),360);
    circstats_slow(i,2) = circ_rad2ang(circ_std(legPhaseDiffs_rads_s(:,i)));
    
end

for i=9:16
    
    circstats_median(i,1) = mod(circ_rad2ang(circ_median(legPhaseDiffs_rads_m(:,i-8))),360);
    circstats_std(i,1) = circ_rad2ang(circ_std(legPhaseDiffs_rads_m(:,i-8)));
    
    circstats_med(i-8,1) = mod(circ_rad2ang(circ_median(legPhaseDiffs_rads_m(:,i-8))),360);
    circstats_med(i-8,2) = circ_rad2ang(circ_std(legPhaseDiffs_rads_m(:,i-8)));
    
end

for i=17:24
    
    circstats_median(i,1) = mod(circ_rad2ang(circ_median(legPhaseDiffs_rads_f(:,i-16))),360);
    circstats_std(i,1) = circ_rad2ang(circ_std(legPhaseDiffs_rads_f(:,i-16)));
    
    circstats_fast(i-16,1) = mod(circ_rad2ang(circ_median(legPhaseDiffs_rads_f(:,i-16))),360);
    circstats_fast(i-16,2) = circ_rad2ang(circ_std(legPhaseDiffs_rads_f(:,i-16)));
    
end

circstats_all = horzcat(circstats_median,circstats_std,circstats_speeds);
circstats_headers = {'leg','median','s.d.','speed'};
circstats_legscol = {'L1';'L2';'L3';'L4';'R1';'R2';'R3';'R4';'L1';'L2';'L3';'L4';'R1';'R2';'R3';'R4';'L1';'L2';'L3';'L4';'R1';'R2';'R3';'R4'};

circstats_final_t = horzcat(circstats_legscol,num2cell(circstats_all));
circstats_final = vertcat(circstats_headers,circstats_final_t);

cell2csv([saveDir,'circstats_legphases_bpN_WolfIntactFESNUnC.csv'], circstats_final);

%% CALCULATE OFFSETS FOR EACH SPEED BIN

% L-R offsets
% - segNum0 R-L R1-L1 (5,1)-(1,1)
% - segNum1 R-L R2-L2 (6,1)-(2,1)
% - segNum2 R-L R3-L3 (7,1)-(3,1)
% - segNum3 R-L R4-L4 (8,1)-(4,1)

offsets_LR_s = nan(4,1);
offsets_LR_m = nan(4,1);
offsets_LR_f = nan(4,1);

for i=1:4
   
    offsets_LR_s(i) = mod(circstats_slow(i+4,1) - circstats_slow(i,1),360);
    offsets_LR_m(i) = mod(circstats_med(i+4,1) - circstats_med(i,1),360);
    offsets_LR_f(i) = mod(circstats_fast(i+4,1) - circstats_fast(i,1),360);
    
end

offsets_LR_all = vertcat(offsets_LR_s,offsets_LR_m,offsets_LR_f);
offsets_LR_headers = {'legs','offset','speed'};
offsets_LR_legscol = {'R1-L1';'R2-L2';'R3-L3';'R4-L4';'R1-L1';'R2-L2';'R3-L3';'R4-L4'; 'R1-L1';'R2-L2';'R3-L3';'R4-L4'};
offsets_LR_speeds = [1;1;1;1;2;2;2;2;3;3;3;3];

offsets_LR_final_t = horzcat(offsets_LR_legscol,num2cell(offsets_LR_all),num2cell(offsets_LR_speeds));
offsets_LR_final = vertcat(offsets_LR_headers,offsets_LR_final_t);

cell2csv([saveDir,'offsets_LR_WolfIntactFESNUnC.csv'], offsets_LR_final);


% SEGMENT OFFSETS
% - L segNum1-segNum0 L2-L1 (2,1)-(1,1)
% - L segNum2-segNum1 L3-L2 (3,1)-(2,1)
% - L segNum3-segNum2 L4-L3 (4,1)-(3,1)

% - R segNum1-segNum0 R2-R1 (6,1)-(5,1)
% - R segNum2-segNum1 R3-R2 (7,1)-(6,1)
% - R segnum3-segNum2 R4-R3 (8,1)-(7,1)

offsets_seg_s = nan(3,2);
offsets_seg_m = nan(3,2);
offsets_seg_f = nan(3,2);


for i=1:3
   
   offsets_seg_s(i,1) = mod(circstats_slow(i+1,1) - circstats_slow(i,1),360);
   offsets_seg_m(i,1) = mod(circstats_med(i+1,1) - circstats_med(i,1),360);
   offsets_seg_f(i,1) = mod(circstats_fast(i+1,1) - circstats_fast(i,1),360);
   
   offsets_seg_s(i,2) = mod(circstats_slow(i+5,1) - circstats_slow(i+4,1),360);
   offsets_seg_m(i,2) = mod(circstats_med(i+5,1) - circstats_med(i+4,1),360);
   offsets_seg_f(i,2) = mod(circstats_fast(i+5,1) - circstats_fast(i+4,1),360);
    
end

offsets_seg_all = vertcat(offsets_seg_s,offsets_seg_m,offsets_seg_m);
offsets_seg_headers = {'segNum','L_offset','R_offset','speed'};
offsets_seg_segcol = {'1-0';'2-1';'3-2';'1-0';'2-1';'3-2';'1-0';'2-1';'3-2'};
offsets_seg_speeds = [1;1;1;2;2;2;3;3;3];

offsets_seg_final_t = horzcat(offsets_seg_segcol,num2cell(offsets_seg_all),num2cell(offsets_seg_speeds));
offsets_seg_final = vertcat(offsets_seg_headers,offsets_seg_final_t);

cell2csv([saveDir,'offsets_seg_WolfIntactFESNUnC.csv'], offsets_seg_final);


% TETRAPOD OFFSETS
% Tetrapod 1: L1 + R2 + L3 + R4
% Offsets: 
% R2-L1 (6,1)-(1,1) | (2)-(1)
% L3-R2 (3,1)-(6,1) | (3)-(2)
% R4-L3 (8,1)-(3,1) | (4)-(3)

% Tetrapod 2: R1 + L2 + R3 + L4
% Offsets: 
% L2-R1 (2,1)-(5,1) | (2)-(1)
% R3-L2 (7,1)-(2,1) | (3)-(2)
% L4-R3 (4,1)-(7,1) | (4)-(3)

tet1_medians_s = [circstats_slow(1,1),circstats_slow(6,1),circstats_slow(3,1),circstats_slow(8,1)];
tet1_medians_m = [circstats_med(1,1),circstats_med(6,1),circstats_med(3,1),circstats_med(8,1)];
tet1_medians_f = [circstats_fast(1,1),circstats_fast(6,1),circstats_fast(3,1),circstats_fast(8,1)];

tet2_medians_s = [circstats_slow(5,1),circstats_slow(2,1),circstats_slow(7,1),circstats_slow(4,1)];
tet2_medians_m = [circstats_med(5,1),circstats_med(2,1),circstats_med(7,1),circstats_med(4,1)];
tet2_medians_f = [circstats_fast(5,1),circstats_fast(2,1),circstats_fast(7,1),circstats_fast(4,1)];

offsets_tet1_s = nan(3,1);
offsets_tet1_m = nan(3,1);
offsets_tet1_f = nan(3,1);

offsets_tet2_s = nan(3,1);
offsets_tet2_m = nan(3,1);
offsets_tet2_f = nan(3,1);

for i=1:3
    
   offsets_tet1_s(i,1) = mod(tet1_medians_s(i+1) - tet1_medians_s(i),360);
   offsets_tet1_m(i,1) = mod(tet1_medians_m(i+1) - tet1_medians_m(i),360);
   offsets_tet1_f(i,1) = mod(tet1_medians_f(i+1) - tet1_medians_f(i),360);
   
   offsets_tet2_s(i,1) = mod(tet2_medians_s(i+1) - tet2_medians_s(i),360);
   offsets_tet2_m(i,1) = mod(tet2_medians_m(i+1) - tet2_medians_m(i),360);
   offsets_tet2_f(i,1) = mod(tet2_medians_f(i+1) - tet2_medians_f(i),360);
   
end

offsets_tet_all = vertcat(offsets_tet1_s,offsets_tet1_m,offsets_tet1_f,offsets_tet2_s,offsets_tet2_m,offsets_tet2_f);
offsets_tet_headers = {'tetNum','legs','offset','speed'};
offsets_tet_legscol = {'R2-L1';'L3-R2';'R4-L3';'R2-L1';'L3-R2';'R4-L3';'R2-L1';'L3-R2';'R4-L3';'L2-R1';'R3-L2';'L4-R3';'L2-R1';'R3-L2';'L4-R3';'L2-R1';'R3-L2';'L4-R3'};
offsets_tet_speeds = [1;1;1;2;2;2;3;3;3;1;1;1;2;2;2;3;3;3];
offsets_tet_tetnum = [1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2];

offsets_tet_final_t = horzcat(num2cell(offsets_tet_tetnum),offsets_tet_legscol,num2cell(offsets_tet_all),num2cell(offsets_tet_speeds));
offsets_tet_final = vertcat(offsets_tet_headers,offsets_tet_final_t);

cell2csv([saveDir,'offsets_tet_WolfIntactFESNUnC.csv'], offsets_tet_final);



