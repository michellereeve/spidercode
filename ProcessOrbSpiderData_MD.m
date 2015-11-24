function [] = ProcessOrbSpiderData()
% A script to process spider gait data (initially orb-weaver, intact trials)

% set ColourOrder for plots - need 8 colours, default is 7
SpidColors = [0 0.4470 0.7410; 
    0.8500 0.3250 0.0980; 
    0.9290 0.6940 0.1250; 
    0.4940 0.1840 0.5560; 
    0.4660 0.6740 0.1880; 
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0 0 0.4828]; 
set(groot,'defaultAxesColorOrder',SpidColors);

% user input data if not already loaded in - uigetfile. Output fileName.


if ~exist('fileName','var')
    %Browse for the file
    [kFilename, kPathname] = uigetfile(['.csv'], 'Choose coordinate data text (csv) file');
else
    kFilename = fileName;
    kPathname = pathName;
end
filePrefix = kFilename(16:(end-7));

%Check for previously processed file, to skip smoothing steps if possible
if ~exist([kPathname filePrefix '.mat'],'file')


% set framerate
framerate = 1000;

% load in data file (tempData) - dlmread / csvread
tempData = csvread([kPathname kFilename]);

% column numbers (for intact data):
% 1 = frame number
% 2 = time (sec)
% 3-4 = bodyCOM X,Y
% 5-6 = bodyBack X,Y
% 7-8 = bodyFront X,Y
% 9-10 = L1 X,Y
% 11-12 = L2 X,Y
% 13-14 = L3 X,Y ... so from #9, odd numbers are X, even numbers are Y
% ends at 24 (R4 Y)

nCols = size(tempData,2);
xLegCols = [9:2:nCols];

% % plot raw data - user choose subset from ginput. Decide if this is X vs time data or calc leg length first
% f1=figure;
% plot(tempData(:,xLegCols));
% cutoff = round(ginput(1));
% close(f1);

% % update data to new subset
% % Empty all the elements at the end of the original vector, from the xcoord you clicked on onwards
% tempData((cutoff(1):end),:) = [];

% ALTERNATE STEP TO USER CHOOSING SUBSET
% remove instances of -1 in kinematic data - usually this means an error in digitising step
isNaN = find(tempData == -1);
tempData(isNaN) = NaN;
%chop the end of the dataset off where the NaNs occur (this screws up Hilbert calcs later on if not removed.) This is a very rough
% way of doing it - if there's a NaN in the middle of the dataset it will
% remove loads of useful data. Should change it to check for >2 adjacent
% rows containing NaNs but will do for now (or can check with Monica!)
[NANROW NANCOL] = find(isnan(tempData));
tempData((NANROW(1):end),:) = [];

% set some variables
time = tempData(:,2);
frames = tempData(:,1);
% row number doesn't change with different subsets once chopped so this variable is always the same
nRows = size(tempData,1);

% Pull out kinematic data
[kineData] = PullOutKinematics (tempData,nRows);

clear tempData
nCols = size(kineData,2); %update to new dataset

% set some variables for kineData - XY data
xLegs = [7:2:size(kineData,2)];
yLegs = [8:2:size(kineData,2)];
xBody = [1:2:6];
yBody = [2:2:6];
xCoords = [1:2:size(kineData,2)]; % all X coords
yCoords = [2:2:size(kineData,2)]; % all Y coords

% spline smooth data - output smData - plot smoothed data
% Run with sptol set to zero to get raw velocity of each column (body XY, leg XY) (rawVel)
[~,rawVel] = SplineInterp_wSPAPS(kineData, time, 0, 1);

% calculate raw resultant velocities - first 3 cols = bodyCOM, bodyBack, bodyFront, then legs 
rawTotVel = sqrt(rawVel(:,xCoords).^2 + rawVel(:,yCoords).^2);

% Run with sptol to get smoothed velocities (smVel)

%Filtering section:  Spline filter data in a loop with user feedback
% Set default tolerance
vel_sptol = 0.2; % sptol for smoothing velocities
SSans=inputdlg('Filter settings (velocities)','Input spline tolerance',1,{num2str(vel_sptol)});
vel_sptol=str2num(SSans{1,1}); 

%Create a 'while' loop, allow user to adjust spline tolerances
filtGood = 'N';
while filtGood == 'N'
    
    [smKineData,smVel] = SplineInterp_wSPAPS(kineData, time, vel_sptol, 1);

    % calculate smoothed resultant velocity (smTotVel) from smoothed XY velocities (smVel)
    smTotVel = sqrt(smVel(:,xCoords).^2 + smVel(:,yCoords).^2);

    % plot some legs raw resultant vel and smoothed resultant vel to check sptol
    f1=figure;
    plot(time,rawTotVel(:,6),'r');
    hold on;
    plot(time, smTotVel(:,6),'b');
    xlabel('Time (s)')
    ylabel('Foot velocity (mm/s)')
    title([filePrefix ': Kinematic data: red = raw, blue = filtered | sptol = ' num2str(vel_sptol) ' CLOSE'])
    [~] = SaveFigAsPDF(f1,kPathname,filePrefix,'_FiltVelocities');
    waitfor(f1);
    
    filtans = inputdlg('Is this filter good? (Y/N)','Filter Status',1,{'N'});
    filtGood=filtans{1,1};
    
    if filtGood == 'N'
        SSans=inputdlg('Adjust Filter','Input new vel spline tolerance',1,{num2str(vel_sptol)});
        vel_sptol=str2num(SSans{1,1}); 
    end
    
end %end while filtering loop for user adjustment of filter settings

% FIGURE OUT WHERE I NEED THESE FUNCTIONS
% [bodyData] = PullOutBodyCoords (kineData,nRows);



% detect foot contacts by velocity thresholds - put in function
[smLegTotVel,smLegTotVelCols,newLegCols] = PullOutLegTotVels (smTotVel,nRows);

% label legs
leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};

%Calculate some velocity values for threshold detection of foot contact
mean_legVel = nanmean(abs(smLegTotVel));
std_legVel = nanstd(abs(smLegTotVel));
% create some empty matrices for the for-loop results
% legContacts = zeros(size(smLegTotVel,1),length(leg_labels));
% legVel_InContact = nan(size(smLegTotVel,1),length(leg_labels));

t_velThreshNum = 0.1; % default velocity threshold number
SSans=inputdlg('Velocity threshold settings','Input velThreshNum',1,{num2str(t_velThreshNum)});
t_velThreshNum=str2num(SSans{1,1}); 

%Create a 'while' loop, allow user to adjust velocity threshold
threshGood = 'N';
while threshGood == 'N'
legContacts = zeros(size(smLegTotVel,1),length(leg_labels));
legVel_InContact = nan(size(smLegTotVel,1),length(leg_labels));
velThreshold = mean_legVel - t_velThreshNum.*(std_legVel);

% For each foot, detect foot contact periods
    for i = 1:8
        velBelowThresh = find(smLegTotVel(:,i) < velThreshold(i));
        legContacts(velBelowThresh,i) = 1;
        legVel_InContact(velBelowThresh,i) = smLegTotVel(velBelowThresh,i);
    end
    
    %Debug plot:
    %Plot the stance phases, overlaid onto kinematic data
    f2=figure;
     hold on;
    for i = 1:4 
    subplot(4,1,i)
    hold on;
    col_i = SpidColors(i,:); 
    col_i2 = SpidColors(i+4,:);    
    plot(time, smLegTotVel(:,i),'Color',col_i);  
    plot(time,legVel_InContact(:,i),'v','Color',col_i)   
    plot(time, smLegTotVel(:,i+4),'Color',col_i2);  
    plot(time,legVel_InContact(:,i+4),'v','Color',col_i2)
    legend(['L ' num2str(i)], 'VD', ['R ' num2str(i)],'VD')
    if i == 2
     ylabel('Leg velocity (mm/s)')
    end  
    end
    xlabel('Time (s)')   
    subplot(4,1,1)
    hold on;
    title([filePrefix ': Detected stance phases, using velocity threshold = ' num2str(t_velThreshNum)])
    [~] = SaveFigAsPDF(f2,kPathname,filePrefix,'_VelocityDetection');
    %close(f2);
    %clear f2; 

    %create gait diagram vectors
    gaitDiagramData = legContacts.* repmat([1 2 3 4 5 6 7 8],length(legContacts),1);
    gaitDiagramData(gaitDiagramData==0) = nan;

    % plot gait diagram
    f3=figure;
    plot(time, gaitDiagramData,'.');
    xlabel('Time (s)')
    ylabel('Stance phases')
    %ylim([0.5 4.5])
    title( [filePrefix ': Initial Gait diagram: Foot velocity only.'])
    set(gca,'YTick',[1 2 3 4 5 6 7 8])
    set(gca,'YTickLabel',leg_labels)
    [~] = SaveFigAsPDF(f3,kPathname,filePrefix,'_GaitDiagram1');
    %waitfor(f3);
    
    threshans = inputdlg('Is this number good? (Y/N)','VelThresh Status',1,{'N'});
    threshGood=threshans{1,1};
    
    if threshGood == 'N'
        SSans=inputdlg('Adjust VelThreshNum','Input new number',1,{num2str(t_velThreshNum)});
        t_velThreshNum=str2num(SSans{1,1}); 
    end
    
    
end %end while loop


% calculate angle of motion & rotation matrix (skip for now)
% plot original & rotated data (skip for now)
% save new XY data (skip for now)

% delete 'over' smoothed data for stance detection. 
clear smKineData

% Use raw kinematic data, smooth again for leg length/angle stuff. Not as
% smoothed as before. Spline filter data in a loop with user feedback
 
% Set default tolerance
kin_sptol = 0.1; % default sptol for calculating kinematic (angle/length) data
SSans=inputdlg('Filter settings (kinetics)','Input spline tolerance',1,{num2str(kin_sptol)});
kin_sptol=str2num(SSans{1,1});

%Create a 'while' loop, allow user to adjust spline tolerances
filtGood = 'N';
while filtGood == 'N'
    
    [filtKineData] = SplineInterp_wSPAPS(kineData, time, kin_sptol, 1);

    % plot some legs raw resultant vel and smoothed resultant vel to check sptol
    f4=figure;
    plot(kineData(:,xLegs),kineData(:,yLegs),'r');
    hold on;
    plot(filtKineData(:,xLegs), filtKineData(:,yLegs),'b');
    xlabel('X-coords')
    ylabel('Y-coords')
    title([filePrefix ': Kinematic data: red = raw, blue = filtered | sptol = ' num2str(kin_sptol) ' CLOSE'])
    [~] = SaveFigAsPDF(f4,kPathname,filePrefix,'_FiltKineData');
    waitfor(f4);
    
    filtans = inputdlg('Is this filter good? (Y/N)','Filter Status',1,{'N'});
    filtGood=filtans{1,1};
    
    if filtGood == 'N'
        SSans=inputdlg('Adjust Filter','Input new kin spline tolerance',1,{num2str(kin_sptol)});
        kin_sptol=str2num(SSans{1,1}); 
    end
    
end %end while filtering loop for user adjustment of filter settings

% ROTATION CODE SHOULD GO HERE, BEOFRE ANYTHING ELSE IS CALCULATED

% pull out body data for subtracting angles from COM
[filtBodyData] = PullOutBodyCoords (filtKineData,nRows);

% pull out just leg data for calculating lengths and angles
[filtLegData, filtXNewLegs, filtYNewLegs] = PullOutLegCoords (filtKineData,nRows);

% Calculate legs Relative to body (Leg X/Y - COM X/Y)
legDataRelToCOM = nan(nRows,size(filtLegData,2));
for i=1:8
    legDataRelToCOM(:,filtXNewLegs(i)) = filtLegData(:,filtXNewLegs(i)) - filtBodyData(:,1);
    legDataRelToCOM(:,filtYNewLegs(i)) = filtLegData(:,filtYNewLegs(i)) - filtBodyData(:,2);
end

% calculate leg lengths & angles
[legLengths,legLengthsMeanSub] = CalcLegLengths (legDataRelToCOM,nRows,filtXNewLegs,filtYNewLegs);
[legAngles,legAnglesMeanSub] = CalcLegAngles (legDataRelToCOM,nRows,filtXNewLegs,filtYNewLegs);

% plot leg lengths and angles for sanity check 
f5 = figure();
plot(time,legLengthsMeanSub);
xlabel('Time (s)');
ylabel('Length (mm)');
legend(leg_labels);
title([filePrefix ': Leg Lengths (rel to COM, mean-subtracted)']);
[~] = SaveFigAsPDF(f5,kPathname,filePrefix,'_LegLengths');
close(f5);

f6 = figure();
for i = 1:4 
subplot(4,1,i)
hold on;
col_i = SpidColors(i,:); 
col_i2 = SpidColors(i+4,:);    
plot(time,legAnglesMeanSub(:,i),'Color',col_i);
plot(time,legAnglesMeanSub(:,i+4),'Color',col_i2);
   legend(['L ' num2str(i)],['R ' num2str(i)])
    if i == 2
     ylabel('Leg Angle (deg)')
    end  
 end
    xlabel('Time (s)')   
    subplot(4,1,1)
    hold on;
    title([filePrefix ': Leg Angles (rel to COM, mean-substracted)'])
    [~] = SaveFigAsPDF(f6,kPathname,filePrefix,'_LegAngles');
    %close(f6);
    
% plot leg orbits - angle vs. length, mean-subtracted
f7 = figure();
plot(legLengthsMeanSub,legAnglesMeanSub);
ylabel('Angle (deg)');
xlabel('Length (mm)');
legend(leg_labels);
title([filePrefix ': Leg Orbital Plot (rel to COM, mean-substracted)']);
[~] = SaveFigAsPDF(f7,kPathname,filePrefix,'_LegOrbitalPlot');
close(f7);

% plot leg orbits - angle vs. length, not mean-subtracted
f8 = figure();
plot(legLengths,legAngles);
ylabel('Angle (deg)');
xlabel('Length (mm)');
legend(leg_labels);
title([filePrefix ': Leg Orbital Plot (rel to COM)']);
[~] = SaveFigAsPDF(f8,kPathname,filePrefix,'_LegOrbitalPlot2');
close(f8);

% think I need to plot the leg orbits on separate plots for each leg, do a
% subplot for each leg. Do I want mean-subtracted or not? Does it matter?

% convert legAngles to radians for the polar plot:
[legAnglesRad] = deg2rad (legAngles);

% plot polar plots. Angle needs to be in radians. The act of calculating the
% angles and lengths (already done) is equivalent to the cart2pol function.
f9=(figure);
polar(legAnglesRad,legLengths);
legend(leg_labels,'Location', 'eastoutside');
title([filePrefix ': Leg Polar Plot (rel to COM)']);
print([kPathname,filePrefix,'_LegPolarPlot'],'-dpdf');
%close(f9);


%Save data after smoothing
save([kPathname filePrefix])

else    
    load([kPathname filePrefix])
end

%%

% calculate Hilbert phases (normal and inverted)
for i=1:8
    hilbert_phase(:,i) = angle(hilbert(legAnglesMeanSub(:,i)));
    hilbert_phase_inverted(:,i) = angle(hilbert(legAnglesMeanSub(:,i)).*-1);
end

%These variables do need to be pre-allocated to ensure
% that they are specified to have the correct size/dimensions
eventStart_H1 = cell(8,1);
eventStart_H2 = cell(8,1);
footDown_I = cell(8,1);
footOff_I = cell(8,1);

stridePeriod = cell(8,1);
stancePeriod = cell(8,1);
swingPeriod = cell(8,1);
dutyFactor = cell(8,1);

newGaitDiagram = nan(size(legContacts));
hilbertEvents = cell(8,1);
velocityEvents = cell(8,1);
combinedEvents = cell(8,1);
refPhase = NaN(size(hilbert_phase,1),8);

for i=1:8
    eventStart_H1{i} = find(diff(hilbert_phase(:,i))<(0.25*pi.*-1));
    eventStart_H2{i} = find(diff(hilbert_phase_inverted(:,i))<(0.25*pi*-1));
    % Right legs: hilbert 1 (H1) events correspond to foot off transitions
    % Left legs: hilbert 2 (H2) events correspond to foot off transitions
    
    % sort rows by eventStart_H1 in ascending order of column 1
    if i <= 4 % Left leg
    hilbertSorted = sortrows([eventStart_H1{i};  eventStart_H2{i} ],1);
    refPhase(:,i) = hilbert_phase_inverted(:,i);
    else %Right leg
    hilbertSorted = sortrows([eventStart_H2{i};  eventStart_H1{i} ],1);
    refPhase(:,i) = hilbert_phase(:,i);
    end
    
    hilbertEvents{i} = hilbertSorted;
    
    footDown_I{i} = find(diff(legContacts(:,i))==1);
    footDown_I{i}(:,2) =  ones(length(footDown_I{i}(:,1)),1);
    
    footOff_I{i} = find(diff(legContacts(:,i))==-1);
    footOff_I{i}(:,2) = ones(length(footOff_I{i}(:,1)),1).*-1;
    
    vthreshSorted = sortrows([ footDown_I{i};  footOff_I{i} ],1);
    velocityEvents{i} = vthreshSorted;
    
    if size(hilbertSorted,1)>=2 && size(vthreshSorted,1)>= 2 %If there are enough events for this leg to calculate full stride cycles,  get data:
        [newEventsSorted]= CombineHilbertVelocityEvents(hilbertSorted,vthreshSorted,smLegTotVel(:,i));
        
        combinedEvents{i} = newEventsSorted;
        halfPeriods_combined = diff(newEventsSorted,1);
        
        stancePeriod{i} = (halfPeriods_combined(halfPeriods_combined(:,2)==-2,1))./framerate;
        swingPeriod{i} = (halfPeriods_combined(halfPeriods_combined(:,2)==2,1))./framerate;
        
        stridePeriod{i} = [diff(newEventsSorted(newEventsSorted(:,2)== 1,1))]./framerate;            

        if isempty(stridePeriod{i})
              stridePeriod{i} = [diff(newEventsSorted(newEventsSorted(:,2)== -1,1))]./framerate;
        end
        
        t_end = min([length(stridePeriod{i}) length(stancePeriod{i})]);
        dutyFactor{i} =  stancePeriod{i}(1:t_end)./ stridePeriod{i}(1:t_end) ;
        
        %New gait diagram
        fD_Events = newEventsSorted(newEventsSorted(:,2)== 1,1);
        fO_Events = newEventsSorted(newEventsSorted(:,2)== -1,1);
        
        if fO_Events(1)<fD_Events(1)
            curStanceIndex = [1:1:fO_Events(1)];
            newGaitDiagram(curStanceIndex,i) = 1.*i;
        end
        
        for k=1:length(fD_Events)
            cFD = fD_Events(k);
            cFO = min([fO_Events(fO_Events>cFD)]);
            curStanceIndex = [cFD:1:min([cFO length(newGaitDiagram)])];
            newGaitDiagram(curStanceIndex,i) = 1.*i;
        end
        
    else %If there are not enough events to calculate full stride cycles for this leg,  enter NaNs
        newEventsSorted = NaN;
        combinedEvents{i} = newEventsSorted;
        stancePeriod{i} = NaN;
        swingPeriod{i} = NaN;
        stridePeriod{i} = NaN;
        dutyFactor{i} =  NaN ;
        
    end
   
    
end

%     f2=figure;
%      hold on;
%     for i = 1:4 
%     subplot(4,1,i)
%     hold on;
%     col_i = SpidColors(i,:); 
%     col_i2 = SpidColors(i+4,:);    
%     plot(time, smLegTotVel(:,i),'Color',col_i);  
%     plot(time,legVel_InContact(:,i),'v','Color',col_i)   
%     plot(time, smLegTotVel(:,i+4),'Color',col_i2);  
%     plot(time,legVel_InContact(:,i+4),'v','Color',col_i2)
%     legend(['L ' num2str(i)], 'VD', ['R ' num2str(i)],'VD')
%     if i == 2
%      ylabel('Leg velocity (mm/s)')
%     end  
%     end
%     xlabel('Time (s)')   
%     subplot(4,1,1)
%     hold on;
%     title([filePrefix ': Detected stance phases, using velocity threshold = ' num2str(t_velThreshNum)])
% 
    
%plot hilbert phases
f10 = figure;
 for i = 1:4 
subplot(4,1,i)
hold on;
col_i = SpidColors(i,:); 
col_i2 = SpidColors(i+4,:);    
plot(time,hilbert_phase(:,i),'Color',col_i);
plot(time,hilbert_phase(:,i+4),'Color',col_i2);
    plot(time(eventStart_H1{i}), hilbert_phase(eventStart_H1{i},i),'kX');
    plot(time(eventStart_H1{i+4}), hilbert_phase(eventStart_H1{i+4},i+4),'kX');
plot(time,hilbert_phase_inverted(:,i),':','Color',col_i);
plot(time,hilbert_phase_inverted(:,i+4),':','Color',col_i2);
    plot(time(eventStart_H2{i}), hilbert_phase_inverted(eventStart_H2{i},i),'kX');
    plot(time(eventStart_H2{i+4}), hilbert_phase_inverted(eventStart_H2{i+4},i+4),'kX');
legend(['H1 L ' num2str(i)],['H1 R ' num2str(i)],'ED', 'ED',['H2 L ' num2str(i)], ['H2 R ' num2str(i)],'ED','ED');
if i == 2
 ylabel('Hilbert phase')
end
 end
xlabel('Time (s)')
subplot(4,1,1)
hold on;
title([filePrefix ': Hilbert phase of leg angles'])
[~] = SaveFigAsPDF(f10,kPathname,filePrefix,'_HilbertPhase');
%close(f10);

 
%Plot the new gait diagram
f11= figure;
plot(time, newGaitDiagram,'.');
    xlabel('Time (s)')
    ylabel('Stance phases')
     set(gca,'YTick',[1 2 3 4 5 6 7 8])
    set(gca,'YTickLabel',leg_labels)
title([filePrefix ': New Gait diagram: foot velocity and leg angle detection'] )
[~] = SaveFigAsPDF(f11,kPathname,filePrefix,'_GaitDiagram_Combined');
    

end


function [kineData, kineCoords, newKineCoords] = PullOutKinematics (tempData,nRows)
% PullOutKinematics takes all kinematics (legs & body) from imported CSV
% data from ProAnalyst, and puts it in a matrix for smoothing.

% index for the columns with kinematic data
kineCoords = [3:size(tempData,2)];

% create empty dataset
kineData = nan(nRows,length(kineCoords));
newKineCoords = [1:size(kineData,2)];

% pull out kinematic data
kineData(:,newKineCoords) = tempData(:,kineCoords);

end


function [bodyData] = PullOutBodyCoords (kineData,nRows)
% PullOutBodyCoords takes body coordinates from imported CSV data from
% ProAnalyst, and puts it in a matrix.

% index for body X and Y columns
xPtsBody = [1:2:6];
yPtsBody = [2:2:6];

% create empty dataset & new XYs
bodyData = nan(nRows,length(xPtsBody).*2);
xNewBody = [1:2:size(bodyData,2)];
yNewBody =  [2:2:size(bodyData,2)];

% pull out body XY data
bodyData(:,xNewBody) = kineData(:,xPtsBody);
bodyData(:,yNewBody) = kineData(:,yPtsBody);
end

function [legData, xNewLegs, yNewLegs] = PullOutLegCoords (kineData,nRows)
% PullOutLegCoords takes leg coordinates from imported CSV data from
% ProAnalyst, and puts it in a matrix.

% index for leg X and Y columns
xPtsLegs = [7:2:size(kineData,2)];
yPtsLegs = [8:2:size(kineData,2)];

% create empty dataset & new XYs (what data type is this??)
legData = nan(nRows,length(xPtsLegs).*2);
xNewLegs = [1:2:size(legData,2)];
yNewLegs =  [2:2:size(legData,2)];

% pull out leg XY data
% in data, all rows with column numbers in xNew - fill with stuff in
% tempData column numbers in xPts
legData(:,xNewLegs) = kineData(:,xPtsLegs);
legData(:,yNewLegs) = kineData(:,yPtsLegs);
end

function [smLegTotVel, smLegTotVelCols,newLegCols] = PullOutLegTotVels (smTotVel,nRows)
% PullOutLegTotVels takes just the leg total velocities from smoothed
% velocities calculated by the SPLINE function. Works on Matrices with 11
% rows (1-3 = body, 4-11 = legs)

% index for leg columns
smLegTotVelCols = [4:size(smTotVel,2)];

% create empty dataset & new columns
smLegTotVel = nan(nRows, length(smLegTotVelCols));
newLegCols = [1:size(smLegTotVel,2)];

% pull out data
smLegTotVel(:,newLegCols) = smTotVel(:,smLegTotVelCols);

end


function [figH] = SaveFigAsPDF(figH,defDir,baseFNameString,suffixString)

if ~exist('suffixString','var')
    suffixString = [];
end
set(figH,'units', 'normalized'); set(figH,'Position', [0 0.0364583 1 0.875]);
figFilename = [defDir baseFNameString suffixString '.pdf'];
saveas(figH,figFilename,'pdf');
end


function [legLengths,legLengthsMeanSub] = CalcLegLengths (legDataRelToCOM,nRows,filtXNewLegs,filtYNewLegs)

legLengths = nan(nRows,length(filtXNewLegs));
legLengthsMeanSub = nan(nRows,length(filtXNewLegs));

for i=1:8
    % calculate leg lengths
    legLengths(:,i) = sqrt((legDataRelToCOM(:,filtXNewLegs(i)).^2)+(legDataRelToCOM(:,filtYNewLegs(i)).^2));
    %subtract mean leg lengths
    legLengthsMeanSub(:,i) = legLengths(:,i) - nanmean(legLengths(:,i));
end

end


function [legAngles,legAnglesMeanSub] = CalcLegAngles (legDataRelToCOM,nRows,filtXNewLegs,filtYNewLegs)

legAngles = nan(nRows, length(filtXNewLegs));
legAnglesMeanSub = nan(nRows,length(filtXNewLegs));

for i=1:8 
    % calculate leg angles
    legAngles(:,i) = atan2d(legDataRelToCOM(:,filtYNewLegs(i)),legDataRelToCOM(:,filtXNewLegs(i)));
    %subtract mean leg angle - NEED TO CHANGE THIS TO CIRCULAR NANMEAN?
    legAnglesMeanSub(:,i) = legAngles(:,i) - nanmean(legAngles(:,i));
%     % unwrapped angles... do I need these?
%     legAnglesMeanSubUnwrap = unwrap(legAnglesMeanSub);
    % find values < 0, replace with 360-abs(angle_value)
    
end


end

function [Rad] = deg2rad (Deg)
% Covert an angle in radians to degrees

Rad = (pi/180) * Deg;

end


function [combinedEvents]= CombineHilbertVelocityEvents(hilbertSorted,vthreshSorted,c_footVel)

newEventsSorted = nan(size(hilbertSorted,1),2);

for      j = 1:size(hilbertSorted,1)
    
    %Check to see if first hilbert event is missing
    % from velocity events (e.g., min diff of first event is
    % larger than min diff for 2nd hilbert event lining  up with 1st velocity event
    % If this is the case,  save the hilbert event in newEventsSorted
    
    c_ind = hilbertSorted(j,1);
    c_diff = abs(vthreshSorted(:,1) - c_ind);
    [c_minVal,c_event] = min(c_diff);
    
    if j < (size(hilbertSorted,1) - 1)
        %check to see if next hilbert event lines up better with the
        %indicated velocity event
        next_ind = hilbertSorted(j+1,1);
        next_diff = abs(vthreshSorted(c_event,1) - next_ind);
        
        if next_diff < c_minVal
            %assume the hilbert event was missing from the velocity events
            % Save the hilbert event as a true event
            newEventsSorted(j,1) = hilbertSorted(j,1);
        else
            newEventsSorted(j,:) = vthreshSorted(c_event,:);
        end
        
    else
        newEventsSorted(j,:) = vthreshSorted(c_event,:);
    end
    
    if j>1
        %             % If the indicated current event is identical to the prior
        %             % event, use hilbert event instead
        %             if  newEventsSorted(j,1) == newEventsSorted(j-1,1)
        %                 newEventsSorted(j,1) = hilbertSorted(j,1);
        %                  newEventsSorted(j,2) = newEventsSorted(j-1,2).*-1;
        %             end
        % If the indicated event is identical to the prior event, insert nan
        if  newEventsSorted(j,1) == newEventsSorted(j-1,1)
            newEventsSorted(j,:) =[NaN NaN];
        end
        
    end
end

if newEventsSorted(end,1) == newEventsSorted(end-1,1)
    newEventsSorted(end,:) =[NaN NaN];
end

newEventsSorted(isnan(newEventsSorted(:,1)),:) = [];

%Polarity of events (footdown/footoff)
%wasn't reliably detected through the combined approache
%Especially for messy velocity data.
% So test for polarity at the end
for v = 1:(size(newEventsSorted,1)-1)
    e1 = newEventsSorted(v,1);
    e2 = newEventsSorted(v+1,1);
    cVel(v) = mean(c_footVel(e1:e2));
end
eventPolarity = sign(diff(cVel));

for v = 1:length(eventPolarity)
    if ~isnan(eventPolarity(v))
        newEventsSorted(v,2) = eventPolarity(v);
    elseif v>1
        newEventsSorted(v,2) = newEventsSorted(v-1,2).*-1;
    end
end
newEventsSorted(end-1,2) = newEventsSorted(end-2,2).*-1;
newEventsSorted(end,2) = newEventsSorted(end-1,2).*-1;
combinedEvents = newEventsSorted;
end

