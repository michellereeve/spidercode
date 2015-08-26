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

% set some variables
time = tempData(:,2);
frames = tempData(:,1);

[kineData, kineCoords, newKineCoords] = PullOutKinematics (tempData);
% PullOutKinematics here - smooth these - body/legs together

% WILL NEED TO CHANGE SPLINE SO IT SMOOTHES PULLOUTKINEMATICS OUTPUT
% spline smooth data - output smData - plot smoothed data
% Run with sptol set to zero to get raw velocity of X and Y (rLegVel)
[~,rLegVel] = SplineInterp_wSPAPS(legData, time, 0, 3);
% calculate resultant raw velocity (rTotVel)
rTotVel = sqrt(rLegVel(:,xNewLegs).^2 + rLegVel(:,yNewLegs).^2);

% Run with sptol to get smoothed velocity of X and Y (legVel)
t_sptol = 0.0025;
[smLegData,legVel] = SplineInterp_wSPAPS(legData, time, t_sptol, 3);
% calculate resultant velocity (totVel) from smoothed XY velocities (legVel)
totVel = sqrt(legVel(:,xNewLegs).^2 + legVel(:,yNewLegs).^2);

% plot raw resultant vel and smoothed resultant vel to check sptol
f1=figure;
    plot(time,rTotVel(:,[1 3]),'r');
    hold on;
    plot(time, totVel(:,[1 3]),'b');
    xlabel('Time')
    ylabel('Foot velocity (mm/s)')
    title([filePrefix ': Kinematic data: red= raw, blue = filtered: CLOSE WINDOW TO CONTINUE'])


% f1=figure;
%     plot(legData(:,xNewLegs), legData(:,yNewLegs),'r');
%     hold on;
%     plot(smLegData(:,xNewLegs), smLegData(:,yNewLegs),'b');
%     xlabel('x-coordinates')
%     ylabel('y-coordinates')
%     title([filePrefix ': Kinematic data: red= raw, blue = filtered: CLOSE WINDOW TO CONTINUE'])

% CHANGE INPUT TO PULLOUTKINEMATICS OUTPUT
[bodyData] = PullOutBodyCoords (tempData);

[legData, xNewLegs, yNewLegs] = PullOutLegCoords (tempData);

clear tempData

% plot raw data - user choose subset from ginput. Decide if this is X vs time data or calc leg length first
f1=figure;
plot(legData(:,xNewLegs));

cutoff = round(ginput(1));

close(f1);

% update data to new subset
% Empty all the elements at the end of the original vector, from the xcoord you clicked on onwards
legData((cutoff(1):end),:) = [];
time((cutoff(1):end),:) = [];

% label legs?
leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};

isNaN = find(legData == -1);
legData(isNaN) = NaN;



% calculate angle of motion & rotation matrix (skip for now)
% plot original & rotated data (skip for now)
% save new XY data (skip for now)

% calculate leg lengths & angles - put in function

% plot leg orbits - angle vs. length

% calculate derivative of xVel & yVel to get overall leg velocity -- got
% this already from SPLINE?

% detect foot contacts by velocity thresholds - put in function

% plot gait diagram

end

function [kineData, kineCoords, newKineCoords] = PullOutKinematics (tempData)
% PullOutKinematics takes all kinematics (legs & body) from imported CSV
% data from ProAnalyst, and puts it in a matrix for smoothing.

% index for the columns with kinematic data
kineCoords = [3:size(tempData,2)];

% create empty dataset
kineData = nan(size(tempData,1),length(kineCoords));
newKineCoords = [1:size(kineData,2)];

% pull out kinematic data
kineData(:,newKineCoords) = tempData(:,kineCoords);

% % index for body X and Y columns - this will be the same across trials
% xPtsBody = [3:2:7];
% yPtsBody = [4:2:8];
% 
% % index for leg X and Y columns - this will work for all trial types
% xPtsLegs = [9:2:size(tempData,2)];
% yPtsLegs = [10:2:size(tempData,2)];
% 
% % create empty dataset & new XYs
% kineData = nan(size(tempData,1),((length(xPtsBody).*2)+(length(xPtsLegs).*2));
% xNewBody = [1:2:size(bodyData,2)];
% yNewBody =  [2:2:size(bodyData,2)];
% xNewLegs = [1:2:size(legData,2)];
% yNewLegs =  [2:2:size(legData,2)];
end


function [bodyData] = PullOutBodyCoords (tempData)
% PullOutBodyCoords takes body coordinates from imported CSV data from
% ProAnalyst, and puts it in a matrix.

% index for body X and Y columns
xPtsBody = [3:2:7];
yPtsBody = [4:2:8];

% create empty dataset & new XYs
bodyData = nan(size(tempData,1),length(xPtsBody).*2);
xNewBody = [1:2:size(bodyData,2)];
yNewBody =  [2:2:size(bodyData,2)];

% pull out body XY data
bodyData(:,xNewBody) = tempData(:,xPtsBody);
bodyData(:,yNewBody) = tempData(:,yPtsBody);
end

function [legData, xNewLegs, yNewLegs] = PullOutLegCoords (tempData)
% PullOutLegCoords takes leg coordinates from imported CSV data from
% ProAnalyst, and puts it in a matrix.

% index for leg X and Y columns
xPtsLegs = [9:2:size(tempData,2)];
yPtsLegs = [10:2:size(tempData,2)];

% create empty dataset & new XYs (what data type is this??)
legData = nan(size(tempData,1),length(xPtsLegs).*2);
xNewLegs = [1:2:size(legData,2)];
yNewLegs =  [2:2:size(legData,2)];

% pull out leg XY data
% in data, all rows with column numbers in xNew - fill with stuff in
% tempData column numbers in xPts
legData(:,xNewLegs) = tempData(:,xPtsLegs);
legData(:,yNewLegs) = tempData(:,yPtsLegs);
end