run importBodyPhaseData.m 

nameString = {'L1','L2','L3', 'L4', 'R1', 'R2','R3', 'R4'};
plotOrd = [1 3 5 7 2 4 6 8];

%reset matlab colour order incase re-running scripts
defaultColorOrder = [ 0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',defaultColorOrder);

% % plot Wolf & Aran separately, just for Intact 
% f0=figure;
% for i = 1:8
% subplot(4,2,plotOrd(i))
% cData = deg2rad([bp_data{bp_WolfIntactFlatInd,(17+i)}]);
% cData(isnan(cData))=[];
% rose(cData,80);
% title(['Leg phase diff intact (flat) ' nameString{i} ' - Wolf=bl Aran=or']) ;
% hold on
% cDataN = deg2rad([bp_data{bp_AranIntactInd,(17+i)}]);
% cDataN(isnan(cDataN))=[];
% rose(cDataN,80);
% end
% %[~] = SaveFigAsPDF_MRPlots(f0,'WolfAranLegPhaseDiffs_flatintact');
% close(f0);

% plot Wolf & Aran separately, just for R4ablt 
% f1=figure;
% for i = 1:8
% subplot(4,2,plotOrd(i))
% cData = deg2rad([bp_data{bp_WolfR4abltFlatInd,(17+i)}]);
% cData(isnan(cData))=[];
% rose(cData,80);
% title(['Leg phase diff R4ablt (flat) ' nameString{i} ' - Wolf=bl Aran=or']) ;
% hold on
% cDataN = deg2rad([bp_data{bp_AranR4abltInd,(17+i)}]);
% cDataN(isnan(cDataN))=[];
% rose(cDataN,80);
% end
% %[~] = SaveFigAsPDF_MRPlots(f1,'WolfAranLegPhaseDiffs_flatR4ablt');
% close(f1);

% % plot Wolf & Aran separately, just for L3aR4m 
% f2=figure;
% for i = 1:8
% subplot(4,2,plotOrd(i))
% cData = deg2rad([bp_data{bp_WolfL3aR4mFlatInd,(17+i)}]);
% cData(isnan(cData))=[];
% rose(cData,80);
% title(['Leg phase diff L3aR4m ' nameString{i} ' - Wolf=bl Aran=or']) ;
% hold on
% cDataN = deg2rad([bp_data{bp_AranL3aR4mInd,(17+i)}]);
% cDataN(isnan(cDataN))=[];
% rose(cDataN,80);
% %[~] = SaveFigAsPDF_MRPlots(f2,'WolfAranLegPhaseDiffs_flatL3aR4m');
% close(f2);
% end


%% STEP 2= MANUAL IMPORT OF CSVs CONTAINING GAIT LINE bp_data AS CELL ARRAY

%% STEP 3 = select subsets
% pull out True/False indexes you need by adding various conditions
bp_WolfIntactFlatESNInd = bp_isWolf & bp_isFlat & bp_isIntact & bp_isESN;
bp_WolfR4abltFlatInd = bp_isWolf & bp_isFlat & bp_isR4ablt & bp_isESN;
bp_WolfL3aR4mFlatInd = bp_isWolf & bp_isFlat & bp_isL3aR4m & bp_isESN;

bp_AranIntactInd = bp_isAran & bp_isFlat & bp_isIntact & bp_isESN;
bp_AranR4abltInd = bp_isAran & bp_isFlat & bp_isR4ablt & bp_isESN;
bp_AranL3aR4mInd = bp_isAran & bp_isFlat & bp_isL3aR4m & bp_isESN;


%% PLOT WOLF bp_data ONLY WITH & WITHOUT GAIT LINES
% plot wolf intact flat no eggsac, no gait lines
f3=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad([bp_data{bp_WolfIntactFlatESNInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' - Wolf intact flat ESN']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f3,'WolfLegPhaseDiffs_allspeeds_intact');
%close(f3)

% plot wolf R4ablt flat, no gait lines
f4=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad([bp_data{bp_WolfR4abltFlatInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' - Wolf R4ablt flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f4,'WolfLegPhaseDiffs_allspeeds_R4ablt');
%close(f4)


% plot wolf L3aR4m flat, no gait lines
f5=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad([bp_data{bp_WolfL3aR4mFlatInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' - Wolf L3aR4m flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f5,'WolfLegPhaseDiffs_allspeeds_L3aR4m');
%close(f5)

% 24/03 chagning to just show alt tetra not metachronal
% set line colours so goes black, blue...
roseWolfLineCols = [0    0   0
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330];
set(groot,'defaultAxesColorOrder',roseWolfLineCols);

% plot wolf intact ESN flat with gait lines
f6=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cAltTet = deg2rad([AltTetWolfFlatESN{:,i}]);
rose(cAltTet,1000);
hold on;
cData = deg2rad([bp_data{bp_WolfIntactFlatESNInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' Wolf intact flat ESN']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f6,'WolfLegPhaseDiffs_allspeeds_intact_AltTet');
%close(f6)


% plot Wolf flat R4ablt with gait lines
f7=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cAltTet = deg2rad([AltTetR4ablt{:,i}]);
rose(cAltTet,1000);
hold on;
cData = deg2rad([bp_data{bp_WolfR4abltFlatInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' Wolf R4ablt flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f7,'WolfLegPhaseDiffs_allspeeds_R4ablt_AltTet');
%close(f7)


% plot Wolf flat L3ablt with gait lines
f8=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cAltTet = deg2rad([AltTetL3aR4m{:,i}]);
rose(cAltTet,1000);
hold on;
cData = deg2rad([bp_data{bp_WolfL3aR4mFlatInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' Wolf L3aR4m flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f8,'WolfLegPhaseDiffs_allspeeds_L3aR4m_AltTet');
%close(f8)


% PLOT ARAN bp_data ONLY WITH & WITHOUT GAIT LINES
% set colours so goes orange...
AranColors = [0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0 0 0.4828];
set(groot,'defaultAxesColorOrder',AranColors);
% plot Aran intact flat, no gait lines
f9=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad([bp_data{bp_AranIntactInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' - Aran intact flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f9,'AranLegPhaseDiffs_allspeeds_intact');
%close(f9)

% plot Aran R4ablt flat, no gait lines
f10=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad([bp_data{bp_AranR4abltInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' - Aran R4ablt flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f10,'AranLegPhaseDiffs_allspeeds_R4ablt');
%close(f10)

% plot Aran L3aR4m flat, no gait lines
f11=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cData = deg2rad([bp_data{bp_AranL3aR4mInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' - Aran L3aR4m flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f11,'AranLegPhaseDiffs_allspeeds_L3aR4m');
%close(f11)

% for gait lines - set line colours so goes black, orange...
roseAranLineCols = [0    0   0
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330];
set(groot,'defaultAxesColorOrder',roseAranLineCols);

% plot Aran intact with gait lines
f12=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cAltTet = deg2rad([AltTetAran{:,i}]);
rose(cAltTet,1000);
hold on;
cData = deg2rad([bp_data{bp_AranIntactInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' Aran intact flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f12,'AranLegPhaseDiffs_allspeeds_intact_AltTet');
%close(f12)

% plot Aran R4ablt with gait lines
f13=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cAltTet = deg2rad([AltTetR4abltAran{:,i}]);
rose(cAltTet,1000);
hold on;
cData = deg2rad([bp_data{bp_AranR4abltInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' Aran R4ablt flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f13,'AranLegPhaseDiffs_allspeeds_R4ablt_AltTet');
%close(f13)

% plot ARan L3ablt with gait lines
f14=figure;
for i = 1:8
subplot(4,2,plotOrd(i))
cAltTet = deg2rad([AltTetL3aR4mAran{:,i}]);
rose(cAltTet,1000);
hold on;
cData = deg2rad([bp_data{bp_AranL3aR4mInd,(17+i)}]);
cData(isnan(cData))=[];
rose(cData,80);
title(['Leg phase diff ' nameString{i} ' Aran L3aR4m flat']) ;
end
%[~] = SaveFigAsPDF_MRPlots(f14,'AranLegPhaseDiffs_allspeeds_L3aR4m_AltTet');
%close(f14)
