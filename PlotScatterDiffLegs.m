% after running createPolyFit

% % T/F matrix of legs for wolf
% lsN_isWolfIntactFlatESN_L1 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegL1;
% lsN_isWolfIntactFlatESN_L2 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegL2;
% lsN_isWolfIntactFlatESN_L3 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegL3;
% lsN_isWolfIntactFlatESN_L4 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegL4;
% lsN_isWolfIntactFlatESN_R1 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegR1;
% lsN_isWolfIntactFlatESN_R2 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegR2;
% lsN_isWolfIntactFlatESN_R3 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegR3;
% lsN_isWolfIntactFlatESN_R4 = lsN_isWolf & lsN_isIntact & lsN_isFlat & lsN_isESN & lsN_isLegR4;

% % T/F matrix of legs for aran
% lsN_isAranIntact_L1 = lsN_isAran & lsN_isIntact & lsN_isLegL1;
% lsN_isAranIntact_L2 = lsN_isAran & lsN_isIntact & lsN_isLegL2;
% lsN_isAranIntact_L3 = lsN_isAran & lsN_isIntact & lsN_isLegL3;
% lsN_isAranIntact_L4 = lsN_isAran & lsN_isIntact & lsN_isLegL4;
% lsN_isAranIntact_R1 = lsN_isAran & lsN_isIntact & lsN_isLegR1;
% lsN_isAranIntact_R2 = lsN_isAran & lsN_isIntact & lsN_isLegR2;
% lsN_isAranIntact_R3 = lsN_isAran & lsN_isIntact & lsN_isLegR3;
% lsN_isAranIntact_R4 = lsN_isAran & lsN_isIntact & lsN_isLegR4;
% 
% %pull out all data for individual legs for wolf
% ls_WolfIntactFlatESN_L1 = lsN_data(lsN_isWolfIntactFlatESN_L1,:);
% ls_WolfIntactFlatESN_L2 = lsN_data(lsN_isWolfIntactFlatESN_L2,:);
% ls_WolfIntactFlatESN_L3 = lsN_data(lsN_isWolfIntactFlatESN_L3,:);
% ls_WolfIntactFlatESN_L4 = lsN_data(lsN_isWolfIntactFlatESN_L4,:);
% ls_WolfIntactFlatESN_R1 = lsN_data(lsN_isWolfIntactFlatESN_R1,:);
% ls_WolfIntactFlatESN_R2 = lsN_data(lsN_isWolfIntactFlatESN_R2,:);
% ls_WolfIntactFlatESN_R3 = lsN_data(lsN_isWolfIntactFlatESN_R3,:);
% ls_WolfIntactFlatESN_R4 = lsN_data(lsN_isWolfIntactFlatESN_R4,:);

%  %pull out all data for individual legs for Aran
% ls_AranIntact_L1 = lsN_data(lsN_isAranIntact_L1,:);
% ls_AranIntact_L2 = lsN_data(lsN_isAranIntact_L2,:);
% ls_AranIntact_L3 = lsN_data(lsN_isAranIntact_L3,:);
% ls_AranIntact_L4 = lsN_data(lsN_isAranIntact_L4,:);
% ls_AranIntact_R1 = lsN_data(lsN_isAranIntact_R1,:);
% ls_AranIntact_R2 = lsN_data(lsN_isAranIntact_R2,:);
% ls_AranIntact_R3 = lsN_data(lsN_isAranIntact_R3,:);
% ls_AranIntact_R4 = lsN_data(lsN_isAranIntact_R4,:);
% 
% % pull out speed data for wolf legs
% ls_WolfIntactFlatESN_L1speed = cell2mat(ls_WolfIntactFlatESN_L1(:,9));
% ls_WolfIntactFlatESN_L2speed = cell2mat(ls_WolfIntactFlatESN_L2(:,9));
% ls_WolfIntactFlatESN_L3speed = cell2mat(ls_WolfIntactFlatESN_L3(:,9));
% ls_WolfIntactFlatESN_L4speed = cell2mat(ls_WolfIntactFlatESN_L4(:,9));
% ls_WolfIntactFlatESN_R1speed = cell2mat(ls_WolfIntactFlatESN_R1(:,9));
% ls_WolfIntactFlatESN_R2speed = cell2mat(ls_WolfIntactFlatESN_R2(:,9));
% ls_WolfIntactFlatESN_R3speed = cell2mat(ls_WolfIntactFlatESN_R3(:,9));
% ls_WolfIntactFlatESN_R4speed = cell2mat(ls_WolfIntactFlatESN_R4(:,9));

% % pull out speed data for Aran legs
% ls_AranIntact_L1speed = cell2mat(ls_AranIntact_L1(:,9));
% ls_AranIntact_L2speed = cell2mat(ls_AranIntact_L2(:,9));
% ls_AranIntact_L3speed = cell2mat(ls_AranIntact_L3(:,9));
% ls_AranIntact_L4speed = cell2mat(ls_AranIntact_L4(:,9));
% ls_AranIntact_R1speed = cell2mat(ls_AranIntact_R1(:,9));
% ls_AranIntact_R2speed = cell2mat(ls_AranIntact_R2(:,9));
% ls_AranIntact_R3speed = cell2mat(ls_AranIntact_R3(:,9));
% ls_AranIntact_R4speed = cell2mat(ls_AranIntact_R4(:,9));
% 
% pull out Slip Factor for wolf legs
% ls_WolfIntactFlatESN_L1slipF = cell2mat(ls_WolfIntactFlatESN_L1(:,23));
% ls_WolfIntactFlatESN_L2slipF = cell2mat(ls_WolfIntactFlatESN_L2(:,23));
% ls_WolfIntactFlatESN_L3slipF = cell2mat(ls_WolfIntactFlatESN_L3(:,23));
% ls_WolfIntactFlatESN_L4slipF = cell2mat(ls_WolfIntactFlatESN_L4(:,23));
% ls_WolfIntactFlatESN_R1slipF = cell2mat(ls_WolfIntactFlatESN_R1(:,23));
% ls_WolfIntactFlatESN_R2slipF = cell2mat(ls_WolfIntactFlatESN_R2(:,23));
% ls_WolfIntactFlatESN_R3slipF = cell2mat(ls_WolfIntactFlatESN_R3(:,23));
% ls_WolfIntactFlatESN_R4slipF = cell2mat(ls_WolfIntactFlatESN_R4(:,23));

% pull out slip factor for Aran
ls_AranIntact_L1slipF = cell2mat(ls_AranIntact_L1(:,23));
ls_AranIntact_L2slipF = cell2mat(ls_AranIntact_L2(:,23));
ls_AranIntact_L3slipF = cell2mat(ls_AranIntact_L3(:,23));
ls_AranIntact_L4slipF = cell2mat(ls_AranIntact_L4(:,23));
ls_AranIntact_R1slipF = cell2mat(ls_AranIntact_R1(:,23));
ls_AranIntact_R2slipF = cell2mat(ls_AranIntact_R2(:,23));
ls_AranIntact_R3slipF = cell2mat(ls_AranIntact_R3(:,23));
ls_AranIntact_R4slipF = cell2mat(ls_AranIntact_R4(:,23));

SpidColors = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0 0 0.4828];
set(groot,'defaultAxesColorOrder',SpidColors);
% plot scatter of speed vs slip for Aran
scatter(ls_AranIntact_L1speed,ls_AranIntact_L1slipF,'filled');
hold on
scatter(ls_AranIntact_L2speed,ls_AranIntact_L2slipF,'filled');
scatter(ls_AranIntact_L3speed,ls_AranIntact_L3slipF,'filled');
scatter(ls_AranIntact_L4speed,ls_AranIntact_L4slipF,'filled');
scatter(ls_AranIntact_R1speed,ls_AranIntact_R1slipF,'filled');
scatter(ls_AranIntact_R2speed,ls_AranIntact_R2slipF,'filled');
scatter(ls_AranIntact_R3speed,ls_AranIntact_R3slipF,'filled');
scatter(ls_AranIntact_R4speed,ls_AranIntact_R4slipF,'filled');
legend('L1','L2','L3','L4','R1','R2','R3','R4');
xlabel('StrideAveVel norm');
ylabel('slipFactor');
title('Aran intact (ls data)');

% plot scatter of speed vs slip for wolfies
scatter(ls_WolfIntactFlatESN_L1speed,ls_WolfIntactFlatESN_L1slipF,'filled');
hold on
scatter(ls_WolfIntactFlatESN_L2speed,ls_WolfIntactFlatESN_L2slipF,'filled');
scatter(ls_WolfIntactFlatESN_L3speed,ls_WolfIntactFlatESN_L3slipF,'filled');
scatter(ls_WolfIntactFlatESN_L4speed,ls_WolfIntactFlatESN_L4slipF,'filled');
scatter(ls_WolfIntactFlatESN_R1speed,ls_WolfIntactFlatESN_R1slipF,'filled');
scatter(ls_WolfIntactFlatESN_R2speed,ls_WolfIntactFlatESN_R2slipF,'filled');
scatter(ls_WolfIntactFlatESN_R3speed,ls_WolfIntactFlatESN_R3slipF,'filled');
scatter(ls_WolfIntactFlatESN_R4speed,ls_WolfIntactFlatESN_R4slipF,'filled');
legend('L1','L2','L3','L4','R1','R2','R3','R4');
xlabel('StrideAveVel norm');
ylabel('slipFactor');
title('Wolf intact flat ESN (ls data)');

% pull out stance Y excur for Aran
ls_AranIntact_L1stanceYExc = cell2mat(ls_AranIntact_L1(:,21));
ls_AranIntact_L2stanceYExc = cell2mat(ls_AranIntact_L2(:,21));
ls_AranIntact_L3stanceYExc = cell2mat(ls_AranIntact_L3(:,21));
ls_AranIntact_L4stanceYExc = cell2mat(ls_AranIntact_L4(:,21));
ls_AranIntact_R1stanceYExc = cell2mat(ls_AranIntact_R1(:,21));
ls_AranIntact_R2stanceYExc = cell2mat(ls_AranIntact_R2(:,21));
ls_AranIntact_R3stanceYExc = cell2mat(ls_AranIntact_R3(:,21));
ls_AranIntact_R4stanceYExc = cell2mat(ls_AranIntact_R4(:,21));

% plot stance Y excur and speed for aran
scatter(ls_AranIntact_L1speed,ls_AranIntact_L1stanceYExc,'filled');
hold on
scatter(ls_AranIntact_L2speed,ls_AranIntact_L2stanceYExc,'filled');
scatter(ls_AranIntact_L3speed,ls_AranIntact_L3stanceYExc,'filled');
scatter(ls_AranIntact_L4speed,ls_AranIntact_L4stanceYExc,'filled');
scatter(ls_AranIntact_R1speed,ls_AranIntact_R1stanceYExc,'filled');
scatter(ls_AranIntact_R2speed,ls_AranIntact_R2stanceYExc,'filled');
scatter(ls_AranIntact_R3speed,ls_AranIntact_R3stanceYExc,'filled');
scatter(ls_AranIntact_R4speed,ls_AranIntact_R4stanceYExc,'filled');
legend('L1','L2','L3','L4','R1','R2','R3','R4');
xlabel('StrideAveVel norm');
ylabel('stance Y excur');
title('Aran intact (ls data)');

% pullt out swing Y excur for aran
ls_AranIntact_L1swingYExc = cell2mat(ls_AranIntact_L1(:,22));
ls_AranIntact_L2swingYExc = cell2mat(ls_AranIntact_L2(:,22));
ls_AranIntact_L3swingYExc = cell2mat(ls_AranIntact_L3(:,22));
ls_AranIntact_L4swingYExc = cell2mat(ls_AranIntact_L4(:,22));
ls_AranIntact_R1swingYExc = cell2mat(ls_AranIntact_R1(:,22));
ls_AranIntact_R2swingYExc = cell2mat(ls_AranIntact_R2(:,22));
ls_AranIntact_R3swingYExc = cell2mat(ls_AranIntact_R3(:,22));
ls_AranIntact_R4swingYExc = cell2mat(ls_AranIntact_R4(:,22));

scatter(ls_AranIntact_L1speed,ls_AranIntact_L1swingYExc,'filled');
hold on
scatter(ls_AranIntact_L2speed,ls_AranIntact_L2swingYExc,'filled');
scatter(ls_AranIntact_L3speed,ls_AranIntact_L3swingYExc,'filled');
scatter(ls_AranIntact_L4speed,ls_AranIntact_L4swingYExc,'filled');
scatter(ls_AranIntact_R1speed,ls_AranIntact_R1swingYExc,'filled');
scatter(ls_AranIntact_R2speed,ls_AranIntact_R2swingYExc,'filled');
scatter(ls_AranIntact_R3speed,ls_AranIntact_R3swingYExc,'filled');
scatter(ls_AranIntact_R4speed,ls_AranIntact_R4swingYExc,'filled');
legend('L1','L2','L3','L4','R1','R2','R3','R4');
xlabel('StrideAveVel norm');
ylabel('swing Y excur');
title('Aran intact (ls data)');


