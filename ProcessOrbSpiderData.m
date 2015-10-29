function [] = ProcessOrbSpiderData()
% A script to process spider gait data (initially orb-weaver, intact trials)

% user input data if not already loaded in - uigetfile. Output fileName.


if ~exist('fileName','var')
    %Browse for the file
    [kFilename, kPathname] = uigetfile(['.csv'], 'Choose coordinate data text (csv) file');
else
    kFilename = fileName;
    kPathname = pathName;
end

filePrefix = kFilename(16:(end-7));

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
% remove instances of -1 - usually this means an error in digitising step
isNaN = find(tempData == -1);
tempData(isNaN) = NaN;

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
vel_sptol = 0.3; % sptol for smoothing velocities
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
legContacts = zeros(size(smLegTotVel,1),length(leg_labels));
legVel_InContact = nan(size(smLegTotVel,1),length(leg_labels));

t_velThreshNum = 0.6; % default velocity threshold number
SSans=inputdlg('Velocity threshold settings','Input velThreshNum',1,{num2str(t_velThreshNum)});
t_velThreshNum=str2num(SSans{1,1}); 

%Create a 'while' loop, allow user to adjust velocity threshold
threshGood = 'N';
while threshGood == 'N'

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
    plot(time, smLegTotVel);
    hold on;
    plot(time,legVel_InContact,'v')
    xlabel('Time (s)')
    ylabel('Leg velocity (mm/s)')
    title([filePrefix ': Detected stance phases, using velocity threshold = ' num2str(t_velThreshNum)])
    [~] = SaveFigAsPDF(f2,kPathname,filePrefix,'_VelocityDetection');
    close(f2);

    %create gait diagram vectors
    gaitDiagramData = legContacts.* repmat([1 2 3 4 5 6 7 8],length(legContacts),1);
    gaitDiagramData(gaitDiagramData==0) = nan;

    % plot gait diagram
    f3=figure;
    plot(time, gaitDiagramData,'.');
    xlabel('Time (s)')
    ylabel('Stance phases')
    %ylim([0.5 4.5])
    title( [filePrefix ': Initial Gait diagram: Foot velocity only. CLOSE'])
    set(gca,'YTick',[1 2 3 4 5 6 7 8])
    set(gca,'YTickLabel',leg_labels)
    [~] = SaveFigAsPDF(f3,kPathname,filePrefix,'_GaitDiagram1');
    waitfor(f3);
    
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
kin_sptol = 0.4; % default sptol for calculating kinematic (angle/length) data
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

% pull out body data for subtracting angles from COM
[filtBodyData] = PullOutBodyCoords (kineData,nRows);

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
plot(time,legAnglesMeanSub);
xlabel('Time (s)');
ylabel('Angle (deg)');
legend(leg_labels);
title([filePrefix ': Leg Angles (rel to COM, mean-substracted)']);
[~] = SaveFigAsPDF(f6,kPathname,filePrefix,'_LegAngles');
close(f6);

% plot leg orbits - angle vs. length, mean-subtracted
f7 = figure();
plot(legLengthsMeanSub,legAnglesMeanSub);
xlabel('Angle (deg)');
ylabel('Length (mm)');
legend(leg_labels);
title([filePrefix ': Leg Orbital Plot (rel to COM, mean-substracted)']);
[~] = SaveFigAsPDF(f7,kPathname,filePrefix,'_LegOrbitalPlot');
close(f7);

% plot leg orbits - angle vs. length, not mean-subtracted
f8 = figure();
plot(legLengths,legAngles);
xlabel('Angle (deg)');
ylabel('Length (mm)');
legend(leg_labels);
title([filePrefix ': Leg Orbital Plot (rel to COM)']);
[~] = SaveFigAsPDF(f8,kPathname,filePrefix,'_LegOrbitalPlot2');
close(f8);

% think I need to plot the leg orbits on separate plots for each leg, do a
% subplot for each leg. Do I want mean-subtracted or not? Does it matter?

% polar plots
% convert cartesian coords to polar using X and Y. Rho = radius (in theory,
% length) and Theta = angle
for i=1:8
    [theta(:,i),rho(:,i)] = cart2pol(legDataRelToCOM(:,filtXNewLegs(i)),legDataRelToCOM(:,filtYNewLegs(i)));
end

%plot polar plots
for i=1:7
    f9 = polar(theta(:,i),rho(:,i));
    hold on
end
polar(theta(:,8),rho(:,8),'k');
legend(leg_labels,'Location', 'eastoutside');
title([filePrefix ': Leg Polar Plot (rel to COM)']);
print([kPathname,filePrefix,'_LegPolarPlot'],'-dpdf');
%save manually if you want it!

% % alternate way of doing polar plots - should give same result as above but
% % doesn't 
% for i=1:7
%     f9 = polar(legAngles(:,i),legLengths(:,i));
%     hold on
% end
% polar(legAngles(:,8),legLengths(:,8),'k');
% legend(leg_labels,'Location', 'eastoutside');
% title([filePrefix ': Leg Polar Plot (rel to COM)']);
% print([kPathname,filePrefix,'_LegPolarPlot'],'-dpdf');

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
