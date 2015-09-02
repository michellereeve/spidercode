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
[~,rawVel] = SplineInterp_wSPAPS(kineData, time, 0, 3);

% calculate raw resultant velocities - first 3 cols = bodyCOM, bodyBack, bodyFront, then legs 
rawTotVel = sqrt(rawVel(:,xCoords).^2 + rawVel(:,yCoords).^2);

% Run with sptol to get smoothed velocities (smVel)

%Filtering section:  Spline filter data in a loop with user feedback
% Set default tolerance
t_sptol = 0.3;
SSans=inputdlg('Filter settings','Input spline tolerance',1,{num2str(t_sptol)});
t_sptol=str2num(SSans{1,1}); 

%Create a 'while' loop, allow user to adjust spline tolerances
filtGood = 'N';
while filtGood == 'N'
    
    [smKineData,smVel] = SplineInterp_wSPAPS(kineData, time, t_sptol, 3);

    % calculate smoothed resultant velocity (smTotVel) from smoothed XY velocities (smVel)
    smTotVel = sqrt(smVel(:,xCoords).^2 + smVel(:,yCoords).^2);

    % plot some legs raw resultant vel and smoothed resultant vel to check sptol
    f1=figure;
        plot(time,rawTotVel(:,4),'r');
        hold on;
        plot(time, smTotVel(:,4),'b');
        xlabel('Time (s)')
        ylabel('Foot velocity (mm/s)')
        title([filePrefix ': Kinematic data: red = raw, blue = filtered: CLOSE WINDOW TO CONTINUE'])
    waitfor(f1);
    
    filtans = inputdlg('Is this filter good? (Y/N)','Filter Status',1,{'N'});
    filtGood=filtans{1,1};
    
    if filtGood == 'N'
        SSans=inputdlg('Adjust Filter','Input new spline tolerance',1,{num2str(t_sptol)});
        t_sptol=str2num(SSans{1,1}); 
    end
    
end %end while filtering loop for user adjustment of filter settings

% FIGURE OUT WHERE I NEED THESE FUNCTIONS
% [bodyData] = PullOutBodyCoords (kineData,nRows);
% [smLegData, smXNewLegs, smYNewLegs] = PullOutLegCoords (smKineData,nRows);


% detect foot contacts by velocity thresholds - put in function
[smLegTotVel,smLegTotVelCols,newLegCols] = PullOutLegTotVels (smTotVel,nRows);

% label legs
leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};

%Calculate some velocity values for threshold detection of foot contact
mean_legVel = nanmean(abs(smLegTotVel));
std_legVel = nanstd(abs(smLegTotVel));
velThreshold = mean_legVel - 0.7.*(std_legVel);

% create some empty matrices for the for-loop results
legContacts = zeros(size(smLegTotVel,1),length(leg_labels));
legVel_InContact = nan(size(smLegTotVel,1),length(leg_labels));

% For each foot, detect foot contact periods
for i = 1:8
    velBelowThresh = find(smLegTotVel(:,i) < velThreshold(i));
    %posBelowThresh = find(rotatedData(:,yPts(allFeet_I(i))) < posThreshold(i));
    %footContacts(intersect(velBelowThresh,posBelowThresh),i) = 1;
    legContacts(velBelowThresh,i) = 1;
    legVel_InContact(velBelowThresh,i) = smLegTotVel(velBelowThresh,i);
end

%Debug plot:
%Plot the stance phases, overlaid onto kinematic data
f3=figure;
plot(time, smLegTotVel);
hold on;
plot(time,legVel_InContact,'v')
xlabel('time (s)')
ylabel('leg velocity (mm/s)')
title([filePrefix ': Detected stance phases, using velocity threshold'])

% plot gait diagram





% calculate angle of motion & rotation matrix (skip for now)
% plot original & rotated data (skip for now)
% save new XY data (skip for now)

% delete 'over' smoothed data for stance detection. Use raw kinematic data, smooth again for leg length/angle stuff. Not as smoothed as before.

% calculate leg lengths & angles - put in function

% plot leg orbits - angle vs. length





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

function [smLegTotVel, smLegTotVelCols,newLegCols] = PullOutLegTotVels (smTotVel,nRows);
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