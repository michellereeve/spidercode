function [out_bp_compiledArray,body_phase_Headers, out_LS_compiledArray, legs_strides_Header, filePrefix, combinedEvents] = ProcessSpiderData_OrbW(fileName,pathName)
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

% user browse to filename if not already provided.

if ~exist('fileName','var')
    %Browse for the file
    [kFilename, kPathname] = uigetfile(['.csv'], 'Choose coordinate data text (csv) file');
else
    kFilename = fileName;
    kPathname = pathName;
end
filePrefix = kFilename(1:(end-13));

%Check for previously processed file, to skip smoothing steps if possible
if ~exist([kPathname filePrefix '.mat'],'file')
    
    % set framerate
    framerate = 1000;
    
    % load in data file (tempData) - dlmread / csvread
    tempData = csvread([kPathname kFilename]);
    
    % Columns (for intact data, without eggsac):
    % 1 = frame number
    % 2 = time (sec)
    % 3-4 = bodyCOM X,Y
    % 5-6 = bodyBack X,Y
    % 7-8 = bodyFront X,Y
    % 9-10 = L1 X,Y
    % 11-12 = L2 X,Y
    % 13-14 = L3 X,Y ... so from #9, odd numbers are X, even numbers are Y
    % ends at 24 (R4 Y)
    
    
    % Columns (for intact data, with eggsac):
    % 1 Frame
    % 2  Time
    % 3-4 bodyCOM X/Y
    % 5-6  bodyBack X/Y
    % 7-8 bodyFront X/Y
    % 9-10
    % eggsac X/Y
    % 11-12 L1 X/Y
    % 13-14 L2 X/Y
    % 15-16 L3 X/Y
    % 17-18   L4 X/Y
    % 19-20 R1 X/Y
    % 21-22 R2 X/Y
    % 23-24 R3 X/Y
    % 25-26 R4 X/Y
    
    
    nCols = size(tempData,2);
    
    %The section below checks for NaN values at the beginning and end
    % of the file and deletes the corresponding rows.
    % It avoids removing NaNs that appear in the middle of the file
    % because these can be interpolated in the filtering steps
    isNaN = find(tempData == -1);
    tempData(isNaN) = NaN; %#ok<*FNDSB>
    
    isNanTempData = isnan(tempData);   
    chckNan = diff(isnan(tempData));    
    isNaNAtStart = isNanTempData(1,:);
    
    nanRow = [];
    for k=1:nCols
        if isNaNAtStart(k)==1
        cRow = find(chckNan(:,k)==-1,1);
        if ~isempty(cRow)
            nanRow = [nanRow; cRow]; %#ok<*AGROW>
        end
        else
            nanRow = [nanRow; 0];
        end
    end
    
    if max(nanRow)~=0
    tempData(1:max(nanRow),:) = [];
    end
    
    isNanTempData = isnan(tempData);   
    chckNan = diff(isnan(tempData));    
    isNaNAtEnd = isNanTempData(end,:);
    
    nanRow = [];
    for k=1:nCols
        if isNaNAtEnd(k)==1
        cRow = find(chckNan(:,k)==1,1,'last');
        if ~isempty(cRow)
            nanRow = [nanRow; cRow];
        end
        else
             nanRow = [nanRow; size(tempData,1)-1];
        end
    end
    
    tempData((min(nanRow)+1):end,:) = [];
    
    % set some variables
    time = tempData(:,2);
    frames = tempData(:,1); %#ok<NASGU>
    % row number doesn't change with different subsets once chopped so this variable is always the same
    nRows = size(tempData,1);
    
    % Pull out kinematic data
    kineCoords = [3:size(tempData,2)];
    kineData = nan(nRows,length(kineCoords));
    newKineCoords = [1:size(kineData,2)];
    kineData(:,newKineCoords) = tempData(:,kineCoords);
    
    clear tempData
    nCols = size(kineData,2); %update to new dataset
    nPts = fix(nCols./2);
    
    %This section below will need to be updated to distinguish 
    % 'intact' and 'ablation' trials
    
%     if strcmp(kFilename(1:2),'03')==1 || strcmp(kFilename(1:2),'04')==1
        % set some variables for kineData - XY data
        xLegsIndex = [7:2:size(kineData,2)];
        yLegsIndex = [8:2:size(kineData,2)];
        xBodyIndex = [1:2:5];
        yBodyIndex = [2:2:6];
%     else
%         % set some variables for kineData - XY data
%         xLegsIndex = [9:2:size(kineData,2)];
%         yLegsIndex = [10:2:size(kineData,2)];
%         xBodyIndex = [1:2:7]; %Additional body point for egg sac
%         yBodyIndex = [2:2:8];
%     end
    
    xCoordsIndex = [1:2:size(kineData,2)]; % all X coords
    yCoordsIndex = [2:2:size(kineData,2)]; % all Y coords
    
    % Reference leg = L2 (Cyclical leg but not adjacent to a later ablated leg
    % or part of the same segment pair)
    refLeg = 2;
    
    % ROTATION CODE
    startIdx = find(sum(isnan(kineData),2) == 0, 1,'first' );
    endIdx = find(sum(isnan(kineData),2) == 0, 1,'last' );
    originPt = kineData(startIdx,1:2); %XY coords of bodyCOM
    kineData = kineData - repmat(originPt,size(kineData,1),nPts);
    
    % Calculate the angle of motion and a rotation matrix
    deltaXtravelled = mean(kineData(endIdx,xBodyIndex)- kineData(startIdx,xBodyIndex));
    deltaYtravelled = mean(kineData(endIdx,yBodyIndex)- kineData(startIdx,yBodyIndex));
    rotation_angle = (atan2(deltaYtravelled,deltaXtravelled)).*-1;
    
    rotation_angle_deg = rad2deg(rotation_angle); %#ok<NASGU>
    rotation_matrix =  [  cos(rotation_angle), -sin(rotation_angle);
        sin(rotation_angle),  cos(rotation_angle)];
    
    rotatedKineData = nan(size(kineData,1),size(kineData,2));
    
    %Rotate the data so that the primary direction of motion is fore-aft
    for k=1:nPts
        c_point= [xCoordsIndex(k)  yCoordsIndex(k)];
        c_vector = kineData(:,c_point)';
        newVector = rotation_matrix*c_vector;
        newVector = newVector';
        rotatedKineData(:,c_point) = newVector;
    end
    
    if  deltaXtravelled <0
        rotatedKineData(:,yPts) = rotatedKineData(:,yPts).*-1;
    end
    
    %Read just origin vertical minimum
    yMin = min(min(rotatedKineData(:,yCoordsIndex)));
    rotatedKineData(:,yCoordsIndex) = rotatedKineData(:,yCoordsIndex) - repmat(yMin,length(rotatedKineData(:,yCoordsIndex)),nPts);
    
    %Plot the original and rotated data for 'reality check'
    f1=figure;
    plot(kineData(:,xCoordsIndex), kineData(:,yCoordsIndex),'r');
    hold on;
    plot(kineData(1,xCoordsIndex), kineData(1,yCoordsIndex),'ro');
    plot(rotatedKineData(:,xCoordsIndex), rotatedKineData(:,yCoordsIndex),'b');
    plot(rotatedKineData(1,xCoordsIndex), rotatedKineData(1,yCoordsIndex),'bo');
    xlabel('x-coordinates')
    ylabel('y-coordinates')
    title([filePrefix ': Rotated kinematic data: red= raw, blue, rotated'])
    [~] = SaveFigAsPDF(f1,kPathname,filePrefix,'_RotatedData');
    close(f1);
    
    % spline smooth data - output smData - plot smoothed data
    % Run with sptol set to zero to get raw velocity of each column (body XY, leg XY) (rawVel)
    [~,rawVel] = SplineInterp_wSPAPS(rotatedKineData, time, 0, 1);
    
    % calculate raw resultant velocities - first 3 cols = bodyCOM, bodyBack, bodyFront, then legs
    rawTotVel = sqrt(rawVel(:,xCoordsIndex).^2 + rawVel(:,yCoordsIndex).^2);
    
    % Run with sptol to get smoothed velocities (smVel)
    % Filtering section:  Spline filter data in a loop with user feedback
    % Set default tolerance
    vel_sptol = 0.001; % sptol for smoothing velocities
    SSans=inputdlg('Filter settings (velocities)','Input spline tolerance',1,{num2str(vel_sptol)});
    vel_sptol=str2num(SSans{1,1});
    
    %Create a 'while' loop, allow user to adjust spline tolerances
    filtGood = 'N';
    while filtGood == 'N'
        
        [~,smVel] = SplineInterp_wSPAPS(rotatedKineData, time, vel_sptol, 1);
        
        % calculate smoothed resultant velocity (smTotVel) from smoothed XY velocities (smVel)
        smTotVel = sqrt(smVel(:,xCoordsIndex).^2 + smVel(:,yCoordsIndex).^2);
        
        % plot some legs raw resultant vel and smoothed resultant vel to check sptol
        f2=figure;
        plot(time,rawTotVel(:,6),'r');
        hold on;
        plot(time, smTotVel(:,6),'b');
        xlabel('Time (s)')
        ylabel('Foot velocity (mm/s)')
        title([filePrefix ': Kinematic data: red = raw, blue = filtered | sptol = ' num2str(vel_sptol) ' CLOSE'])
        [~] = SaveFigAsPDF(f2,kPathname,filePrefix,'_FiltVelocities');
        waitfor(f2);
        
        filtans = inputdlg('Is this filter good? (Y/N)','Filter Status',1,{'N'});
        filtGood=filtans{1,1};
        
        if filtGood == 'N'
            SSans=inputdlg('Adjust Filter','Input new vel spline tolerance',1,{num2str(vel_sptol)});
            vel_sptol=str2num(SSans{1,1});
        end
        
    end %end while filtering loop for user adjustment of filter settings
    
    % detect foot contacts by velocity thresholds
    smLegTotVel = smTotVel(:,(yLegsIndex/2));
    
    % label legs
    leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};
    leg_labels_anat = {'R4','L4','R3','L3','R2','L2','R1','L1'};
    
    %Calculate some velocity values for threshold detection of foot contact
    mean_legVel = nanmean(abs(smLegTotVel));
    std_legVel = nanstd(abs(smLegTotVel));
    % create some empty matrices for the for-loop results
    % legContacts = zeros(size(smLegTotVel,1),length(leg_labels));
    % legVel_InContact = nan(size(smLegTotVel,1),length(leg_labels));
    
    t_velThreshNum = 0.5; % default velocity threshold number
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
        f3=figure;
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
        [~] = SaveFigAsPDF(f3,kPathname,filePrefix,'_VelocityDetection');
        
        %create gait diagram vectors
        gaitDiagramData = legContacts.* repmat([8 6 4 2 7 5 3 1],length(legContacts),1);
        gaitDiagramData(gaitDiagramData==0) = nan;
        
        % plot gait diagram - NEW - in anatomical pairs as before
        f4=figure;
        hold on;
        for i=1:4
            col_i = SpidColors(i,:);
            col_i2 = SpidColors(i+4,:);
            plot(time, gaitDiagramData(:,i),'.','Color',col_i);
            plot(time, gaitDiagramData(:,i+4),'.','Color',col_i2);
        end
        xlabel('Time (s)')
        ylabel('Stance phases')
        title( [filePrefix ': Initial Gait diagram: Foot velocity only.'])
        set(gca,'YTick',[1 2 3 4 5 6 7 8])
        set(gca,'YTickLabel',leg_labels_anat)
        [~] = SaveFigAsPDF(f4,kPathname,filePrefix,'_GaitDiagram1');
        
        threshans = inputdlg('Is this number good? (Y/N)','VelThresh Status',1,{'N'});
        threshGood=threshans{1,1};
        
        if threshGood == 'N'
            SSans=inputdlg('Adjust VelThreshNum','Input new number',1,{num2str(t_velThreshNum)});
            t_velThreshNum=str2num(SSans{1,1});
        end
        
        
    end %end while loop
    
    % delete 'over' smoothed data for stance detection.
    clear smKineData
    
    % Use raw kinematic data, smooth again for leg length/angle stuff. Not as
    % smoothed as before. Spline filter data in a loop with user feedback
    
    % Set default tolerance to same as that of velocity. Can change for a
    % separate smoothing step if I feel I need it in later files/analysis.
    kin_sptol = vel_sptol; % default sptol for calculating kinematic (angle/length) data
    
    [filtRotKineData] = SplineInterp_wSPAPS(rotatedKineData, time, kin_sptol, 1);
    
    % plot some legs raw resultant vel and smoothed resultant vel to check sptol
    f5=figure;
    plot(rotatedKineData(:,xLegsIndex),rotatedKineData(:,yLegsIndex),'r');
    hold on;
    plot(filtRotKineData(:,xLegsIndex), filtRotKineData(:,yLegsIndex),'b');
    xlabel('X-coords')
    ylabel('Y-coords')
    title([filePrefix ': Kinematic data: red = raw, blue = filtered | sptol = ' num2str(kin_sptol)])
    [~] = SaveFigAsPDF(f5,kPathname,filePrefix,'_FiltKineData');
    %close(f5);
    
    % % pull out body data for subtracting angles from COM
    % [filtRBodyData] = PullOutBodyCoords(filtRotKineData,xBodyIndex,yBodyIndex);
    
    % pull out just leg data for calculating lengths and angles
    [filtRLegData, filtXLegCols, filtYLegCols] = PullOutLegCoords(filtRotKineData,xLegsIndex,yLegsIndex);
    
    % Calculate legs Relative to body (Leg X/Y - COM X/Y)
    legDataRelToCOM = nan(nRows,size(filtRLegData,2));
    for i=1:8
        legDataRelToCOM(:,filtXLegCols(i)) = filtRLegData(:,filtXLegCols(i)) - mean(filtRotKineData(:,xBodyIndex),2);
        legDataRelToCOM(:,filtYLegCols(i)) = filtRLegData(:,filtYLegCols(i)) - mean(filtRotKineData(:,yBodyIndex),2);
    end
    
    % % calculate leg lengths & angles
    % [legLengths] = CalcLegLengths(legDataRelToCOM,nRows,filtXLegCols,filtYLegCols);
    % [legAngles] = CalcLegAngles(legDataRelToCOM,nRows,filtXLegCols,filtYLegCols);
    
    aveLegX = nan(1,length(filtXLegCols));
    aveLegY = nan(1,length(filtYLegCols));
    aveLegLength = nan(1,length(filtYLegCols));
    
    % Calculate legs Relative to body (Leg X/Y - COM X/Y)
    legDataRelToCentroid = nan(nRows,size(filtRLegData,2));
    for i=1:8
        aveLegX(i) = mean(legDataRelToCOM(:,filtXLegCols(i)));
        aveLegY(i) = mean(legDataRelToCOM(:,filtYLegCols(i)));
        
        legDataRelToCentroid(:,filtXLegCols(i)) = legDataRelToCOM(:,filtXLegCols(i)) - aveLegX(i);
        legDataRelToCentroid(:,filtYLegCols(i)) = legDataRelToCOM(:,filtYLegCols(i)) - aveLegY(i);
        aveLegLength(i) = sqrt((mean(legDataRelToCOM(:,filtXLegCols(i)))).^2 +(mean(legDataRelToCOM(:,filtYLegCols(i)))).^2);
    end
    
    %     leg_X_RelToCom = legDataRelToCOM(:,filtXLegCols);
%     leg_Y_RelToCom = legDataRelToCOM(:,filtYLegCols);

    leg_X_RelToCom = legDataRelToCentroid(:,filtXLegCols);
    leg_Y_RelToCom = legDataRelToCentroid(:,filtYLegCols);
    
    legDataRelToCentroid(:,filtXLegCols) = legDataRelToCentroid(:,filtXLegCols)+ aveLegX(refLeg);
    legDataRelToCentroid(:,filtYLegCols) = legDataRelToCentroid(:,filtYLegCols)+ aveLegY(refLeg);
    
    % calculate leg lengths & angles
    [legLengthsCentroid,legLengthsMeanSub] = CalcLegLengths(legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols); %#ok<ASGLU>
    [legAnglesCentroid,legAnglesMeanSub] = CalcLegAngles(legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols); %#ok<ASGLU>

    %Save filtered data in  Mat file
    close all; %Close figures to avoid saving them in the mat file
    save([kPathname filePrefix])
    
else    %If data had previously been filtered and saved (.mat file exists)
    kCurrentPath = kPathname; %Ensure that current path is not overwritten
    load([kPathname filePrefix])
    kPathname = kCurrentPath; %Make sure current path is used.
    clear kCurrentPath;
end

% plot foot x y displacements for sanity check
flx = figure;
for i=1:4
    subplot(4,1,i)
    hold on;
    col_i = SpidColors(i,:);
    col_i2 = SpidColors(i+4,:);
    plot(time,leg_X_RelToCom(:,i),'Color',col_i);
    plot(time,leg_X_RelToCom(:,(i+4)),'Color',col_i2);
    legend(['L ' num2str(i)],['R ' num2str(i)])
    if i == 2
        ylabel('Foot x displacement (mm)')
    end
end
xlabel('Time (s)');
subplot(4,1,1)
hold on;
title([filePrefix ': Foot x displacement (rel to COM)']);
[~] = SaveFigAsPDF(flx,kPathname,filePrefix,'_Foot_XDisp');
%close(flx);

% plot foot x y displacements for sanity check
flx = figure;
for i=1:4
    subplot(4,1,i)
    hold on;
    col_i = SpidColors(i,:);
    col_i2 = SpidColors(i+4,:);
    plot(time,leg_Y_RelToCom(:,i),'Color',col_i);
    plot(time,leg_Y_RelToCom(:,(i+4)),'Color',col_i2);
    legend(['L ' num2str(i)],['R ' num2str(i)])
    if i == 2
        ylabel('Foot x displacement (mm)')
    end
end
xlabel('Time (s)');
subplot(4,1,1)
hold on;
title([filePrefix ': Foot y displacement (rel to COM)']);
[~] = SaveFigAsPDF(flx,kPathname,filePrefix,'_Foot_YDisp');
%close(flx);

% % plot leg lengths and angles for sanity check
% f6 = figure();
% for i=1:4
%     subplot(4,1,i)
%     hold on;
%     col_i = SpidColors(i,:);
%     col_i2 = SpidColors(i+4,:);
%     plot(time,legLengthsMeanSub(:,i),'Color',col_i);
%     plot(time,legLengthsMeanSub(:,i+4),'Color',col_i2);
%     legend(['L ' num2str(i)],['R ' num2str(i)])
%     if i == 2
%         ylabel('Leg length (mm)')
%     end
% end
% xlabel('Time (s)');
% subplot(4,1,1)
% hold on;
% title([filePrefix ': Leg Lengths (rel to COM motion & foot centroid)']);
% [~] = SaveFigAsPDF(f6,kPathname,filePrefix,'_LegLengths');
% close(f6);

% f7 = figure();
% for i = 1:4
%     subplot(4,1,i)
%     hold on;
%     col_i = SpidColors(i,:);
%     col_i2 = SpidColors(i+4,:);
%     plot(time,legAnglesMeanSub(:,i),'Color',col_i);
%     plot(time,legAnglesMeanSub(:,i+4),'Color',col_i2);
%     legend(['L ' num2str(i)],['R ' num2str(i)])
%     if i == 2
%         ylabel('Leg Angle (deg)')
%     end
% end
% xlabel('Time (s)')
% subplot(4,1,1)
% hold on;
% title([filePrefix ': Leg Angles (rel to COM motion & foot centroid)'])
% [~] = SaveFigAsPDF(f7,kPathname,filePrefix,'_LegAngles');
% close(f7);

%Plot 'overhead' view plot in cartesian coordinates
% Replaces previous polar plot
f10=(figure); %#ok<NASGU>
plot(legDataRelToCOM(:,filtXLegCols),legDataRelToCOM(:,filtYLegCols));
legend(leg_labels,'Location', 'eastoutside');
title([filePrefix ': Foot trajectory (rel to COM)']);
print([kPathname,filePrefix,'_FootTrajectoryPlot'],'-dpdf');
%close(f10);

% calculate Hilbert phases (normal and inverted)
for i=1:8
    hilbert_phase(:,i) = angle(hilbert(leg_X_RelToCom(:,i)));
    hilbert_phase_inverted(:,i) = angle(hilbert(leg_X_RelToCom(:,i)).*-1);
end

%These variables  need to be pre-allocated to ensure
% that they are specified to have the correct size/dimensions
footOff_H1 = cell(8,1);
footDown_H2 = cell(8,1);
footDown_I = cell(8,1);
footOff_I = cell(8,1);

stridePeriod = cell(8,1);
stancePeriod = cell(8,1);
swingPeriod = cell(8,1);
dutyFactor = cell(8,1);
stance_x_Excur = cell(8,1);
swing_x_Excur = cell(8,1);
stance_y_Excur = cell(8,1);
swing_y_Excur = cell(8,1);
stanceSlipFactor = cell(8,1);

refPhase = NaN(size(hilbert_phase,1),8);

% set up for hilbert event detection while loop for user input
hilEventDet_thresh = 0.06; % default threshold for detecting events from Hilbert
SSans=inputdlg('HilbertEvent threshold settings','Input hilEventDet_thresh',1,{num2str(hilEventDet_thresh)});
hilEventDet_thresh=str2num(SSans{1,1});

%Create a 'while' loop, allow user to adjust hilbert threshold
threshGood = 'N';
while threshGood == 'N'
    %Need to recreate these variables in each iteration of while loop
    %To ensure they are recalculated correctly if the hilbert threshold is adjusted
    newGaitDiagram = nan(size(legContacts));
    hilbertEvents = cell(8,1);
    velocityEvents = cell(8,1);
    combinedEvents = cell(8,1);
    
    for i=1:8
        footOff_H1{i} = find(diff(hilbert_phase(:,i))<(hilEventDet_thresh*pi*-1)); %Foot off events
        footDown_H2{i} = find(diff(hilbert_phase_inverted(:,i))<(hilEventDet_thresh*pi*-1)); %Foot down events

        hilbertSorted = sortrows([footDown_H2{i}; footOff_H1{i}],1);
        refPhase(:,i) = hilbert_phase(:,i);
        
        hilbertEvents{i} = hilbertSorted;
        
        footDown_I{i} = find(diff(legContacts(:,i))==1);
        footDown_I{i}(:,2) =  ones(length(footDown_I{i}(:,1)),1);
        
        footOff_I{i} = find(diff(legContacts(:,i))==-1);
        footOff_I{i}(:,2) = ones(length(footOff_I{i}(:,1)),1).*-1;
        
        vthreshSorted = sortrows([footDown_I{i}; footOff_I{i}],1);
        velocityEvents{i} = vthreshSorted;
        
        if size(hilbertSorted,1)>=2 && size(vthreshSorted,1)>= 2 %If there are enough events for this leg to calculate full stride cycles,  get data:
            [newEventsSorted]= CombineHilbertVelocityEvents(hilbertSorted,vthreshSorted,smLegTotVel(:,i));
            
            combinedEvents{i} = newEventsSorted;
            halfPeriods_combined = diff(newEventsSorted,1);
            
            stancePeriod{i} = (halfPeriods_combined(halfPeriods_combined(:,2)==-2,1))./framerate;
            swingPeriod{i} = (halfPeriods_combined(halfPeriods_combined(:,2)==2,1))./framerate;
            
            % indexing to adjacent '-1' events, which correspond to foot off
            stridePeriod{i} = [diff(newEventsSorted(newEventsSorted(:,2)== -1,1))]./framerate;
            
            t_end = min([length(stridePeriod{i}) length(stancePeriod{i})]);
            dutyFactor{i} =  stancePeriod{i}(1:t_end)./ stridePeriod{i}(1:t_end);
            
            % find values where dutyFactor is <0 - this indicates some
            % problematic stride/stance periods. Replace these strides with NaN
            badValues = find([dutyFactor{i,:}]<0 | [dutyFactor{i,:}]>1);
            dutyFactor{i,1}(badValues) = NaN;%Note that the index badValues comes after the cell index.
            stancePeriod{i,1}(badValues) = NaN;
            swingPeriod{i,1}(badValues) = NaN;
            stridePeriod{i,1}(badValues) = NaN;
            
            %New gait diagram
            fD_Events = newEventsSorted(newEventsSorted(:,2)== 1,1);
            fO_Events = newEventsSorted(newEventsSorted(:,2)== -1,1);
            
            if fO_Events(1)<fD_Events(1)
                curStanceIndex = [1:1:fO_Events(1)];
                newGaitDiagram(curStanceIndex,i) = 1;%.*i;
            end
            
            for k=1:length(fD_Events)
                cFD = fD_Events(k);
                cFO = min([fO_Events(fO_Events>cFD)]); %#ok<*NBRAK>
                curStanceIndex = [cFD:1:min([cFO length(newGaitDiagram)])];
                newGaitDiagram(curStanceIndex,i) = 1;%.*i;
            end
            
%             if i == 8
%                 disp(fO_Events)
%                 disp(fD_Events)
%             end
            
            % calculate leg angle excursions - needs editing to avoid
            % duplicate/extra events, talk to MD
            Sn = min(length(fD_Events),length(fO_Events));
            for S = 1:Sn
                %Calculations need to be adjusted depending on whether swing or stance
                % occur first in the recorded sequence, while minimising 'lost' data
                %If a stance phase starts
                if fO_Events(S,1)>fD_Events(S,1)
                    c_foot_X_stance = leg_X_RelToCom(fD_Events(S,1):fO_Events(S,1),i);
                    c_foot_Y_stance = leg_Y_RelToCom(fD_Events(S,1):fO_Events(S,1),i);
                    %Calculate 'slip factor'
                    c_stanceSlipFactor = sum(~isnan(gaitDiagramData(fD_Events(S,1):fO_Events(S,1),i)))./sum(~isnan(newGaitDiagram(fD_Events(S,1):fO_Events(S,1),i)));
                    if S<= length(fD_Events)-1;
                        c_foot_X_swing = leg_X_RelToCom(fO_Events(S,1):fD_Events(S+1,1),i);
                        c_foot_Y_swing = leg_Y_RelToCom(fO_Events(S,1):fD_Events(S+1,1),i);
                    else
                        c_foot_X_swing = [];
                        c_foot_Y_swing = [];
                    end
                    
                    %if a swing phase starts
                else
                    if S<= length(fO_Events)-1;
                        c_foot_X_stance = leg_X_RelToCom(fD_Events(S,1):fO_Events(S+1,1),i);
                        c_foot_Y_stance = leg_Y_RelToCom(fD_Events(S,1):fO_Events(S+1,1),i);
                        c_stanceSlipFactor = sum(~isnan(gaitDiagramData(fD_Events(S,1):fO_Events(S+1,1),i)))./sum(~isnan(newGaitDiagram(fD_Events(S,1):fO_Events(S+1,1),i)));
                        
                    else
                        c_foot_X_stance = [];
                        c_foot_Y_stance = [];
                    end
                    c_foot_X_swing = leg_X_RelToCom(fO_Events(S,1):fD_Events(S,1),i);
                    c_foot_Y_swing = leg_Y_RelToCom(fO_Events(S,1):fD_Events(S,1),i);
                end
                
                if isempty(c_foot_Y_stance)
                    stance_x_Excur{i,1}(S,1) = NaN;
                    stance_y_Excur{i,1}(S,1)= NaN;
                    stanceSlipFactor{i,1}(S,1) = NaN;
                else
                    stance_x_Excur{i,1}(S,1) = max(c_foot_X_stance)-min(c_foot_X_stance);
                    stance_y_Excur{i,1}(S,1)= max(c_foot_Y_stance)-min(c_foot_Y_stance);
                    stanceSlipFactor{i,1}(S,1) = c_stanceSlipFactor;
                end
                
                if isempty(c_foot_Y_swing)
                    swing_x_Excur{i,1}(S,1) = NaN;
                    swing_y_Excur{i,1}(S,1)= NaN;
                else
                    swing_x_Excur{i,1}(S,1) =  max(c_foot_X_swing)-min(c_foot_X_swing);
                    swing_y_Excur{i,1}(S,1)= max(c_foot_Y_swing)-min(c_foot_Y_swing);
                end
                
            end
            
        else %If there are not enough events to calculate full stride cycles for this leg,  enter NaNs
            newEventsSorted = NaN;
            combinedEvents{i} = newEventsSorted;
            stancePeriod{i} = NaN;
            swingPeriod{i} = NaN;
            stridePeriod{i} = NaN;
            dutyFactor{i} =  NaN ;
            swing_x_Excur{i} =  NaN;
            swing_y_Excur{i} =  NaN;
            stance_x_Excur{i} =  NaN;
            stance_y_Excur{i} =  NaN;
            stanceSlipFactor{i} =  NaN;
        end
    end
    
    %plot reference phases
    f14 = figure;
    for i = 1:4
        subplot(4,1,i)
        hold on;
        col_i = SpidColors(i,:);
        col_i2 = SpidColors(i+4,:);
        plot(time,refPhase(:,i),'Color',col_i);
        plot(time,refPhase(:,i+4),'Color',col_i2);
        % plot events for left legs (H2 = foot off)
        plot(time(footDown_H2{i}), refPhase(footDown_H2{i},i),'ko'); % foot off
        plot(time(footOff_H1{i}), refPhase(footOff_H1{i},i),'kX'); %foot on
        % plot events for right legs (H1 = foot off)
        plot(time(footOff_H1{i+4}), refPhase(footOff_H1{i+4},i+4),'ko'); % foot off
        plot(time(footDown_H2{i+4}), refPhase(footDown_H2{i+4},i+4),'kX'); % foot on
        
        legend(['L' num2str(i)],['R' num2str(i)],'FO', 'FD');
        if i == 2
            ylabel('Hilbert phase')
        end
    end
    xlabel('Time (s)')
    subplot(4,1,1)
    hold on;
    title([filePrefix ': Reference phases (with foot off/down events) | HilbertEvent threshold = ' num2str(hilEventDet_thresh) ' CLOSE'])
    [~] = SaveFigAsPDF(f14,kPathname,filePrefix,'_refPhase');
    waitfor(f14);
    
%     %plot hilbert phases
%     f11 = figure;
%     for i = 1:4
%         subplot(4,1,i)
%         hold on;
%         col_i = SpidColors(i,:);
%         col_i2 = SpidColors(i+4,:);
%         plot(time,hilbert_phase(:,i),'Color',col_i);
%         plot(time,hilbert_phase(:,i+4),'Color',col_i2);
%         plot(time(footOff_H1{i}), hilbert_phase(footOff_H1{i},i),'kX');
%         plot(time(footOff_H1{i+4}), hilbert_phase(footOff_H1{i+4},i+4),'kX');
%         plot(time,hilbert_phase_inverted(:,i),'--','Color',col_i);
%         plot(time,hilbert_phase_inverted(:,i+4),'--','Color',col_i2);
%         plot(time(footDown_H2{i}), hilbert_phase_inverted(footDown_H2{i},i),'kX');
%         plot(time(footDown_H2{i+4}), hilbert_phase_inverted(footDown_H2{i+4},i+4),'kX');
%         legend(['H1 L ' num2str(i)],['H1 R ' num2str(i)],'ED', 'ED',['H2 L ' num2str(i)], ['H2 R ' num2str(i)],'ED','ED');
%         if i == 2
%             ylabel('Hilbert phase')
%         end
%     end
%     xlabel('Time (s)')
%     subplot(4,1,1)
%     hold on;
%     title([filePrefix ': Hilbert phase of leg angles'])
%     [~] = SaveFigAsPDF(f11,kPathname,filePrefix,'_HilbertPhase');
%     close(f11);
    
    %Plot the new gait diagram
    newGaitDiagramData = newGaitDiagram.* repmat([8 6 4 2 7 5 3 1],length(newGaitDiagram),1);
    f12= figure;
    hold on
    for i=1:4
        col_i = SpidColors(i,:);
        col_i2 = SpidColors(i+4,:);
        plot(time, newGaitDiagramData(:,i),'.','Color',col_i);
        plot(time, newGaitDiagramData(:,i+4),'.','Color',col_i2);
    end
    xlabel('Time (s)')
    ylabel('Stance phases')
    set(gca,'YTick',[1 2 3 4 5 6 7 8])
    set(gca,'YTickLabel',leg_labels_anat)
    title([filePrefix ': New Gait diagram: foot velocity and leg angle detection'] )
    [~] = SaveFigAsPDF(f12,kPathname,filePrefix,'_GaitDiagram_Combined');
    
    threshans = inputdlg('Hilbert thresh OK? (Y/N)','hilEventDet_thresh Status',1,{'N'});
    threshGood=threshans{1,1};
    
    if threshGood == 'N'
        SSans=inputdlg('Adjust hilEventDet_thresh','Input new number',1,{num2str(hilEventDet_thresh)});
        hilEventDet_thresh=str2num(SSans{1,1}); %#ok<*ST2NM>
    end
    
end %end while loop

% SAVE OUT STRIDE PARAMETERS AS CSV - LEGS/STRIDES

%Create empty matrix for final leg matrix (all legs compiled)
legs_strides_SaveMatrix = [];

for i=1:8 %For all 8 legs
    c_nRws = max([length(stancePeriod{i}) length(swingPeriod{i}) length(stridePeriod{i})]);
    %Matrix for current leg data
    legs_strides_tempMatrix = nan(c_nRws,11); %create temp matrix full of NaNs
    
    legs_strides_tempMatrix(1:c_nRws,1) = repmat(i,c_nRws,1); %leg number
    legs_strides_tempMatrix(1:c_nRws,2) = [1:c_nRws]; %stride number
    legs_strides_tempMatrix(1:(length(stridePeriod{i})),3) = stridePeriod{i}; %note indexing here in case this variable has less than the max strides
    legs_strides_tempMatrix(1:(length(swingPeriod{i})),4) = swingPeriod{i};
    legs_strides_tempMatrix(1:(length(stancePeriod{i})),5) = stancePeriod{i};
    legs_strides_tempMatrix(1:(length(dutyFactor{i})),6) = dutyFactor{i};
    legs_strides_tempMatrix(1:(length(stance_x_Excur{i})),7) = stance_x_Excur{i};
    legs_strides_tempMatrix(1:(length(swing_x_Excur{i})),8) = swing_x_Excur{i};
    legs_strides_tempMatrix(1:(length(stance_y_Excur{i})),9) = stance_y_Excur{i};
    legs_strides_tempMatrix(1:(length(swing_y_Excur{i})),10) = swing_y_Excur{i};
    legs_strides_tempMatrix(1:(length(stanceSlipFactor{i})),11) = stanceSlipFactor{i};
    
    %concatenate current leg data into the final leg matrix
    legs_strides_SaveMatrix = [legs_strides_SaveMatrix; legs_strides_tempMatrix];
end

%Create a cell array to save as a CSV that includes filename information and headers within the file.
%Filename is repeated in each row so that multiple files can be compiled into a single spreadsheet later.

%Create output header row
legs_strides_Header = {'fileName' 'legNumber', 'strideNumber', 'stridePeriod', 'swingPeriod', ...
    'stancePeriod','dutyFactor','stance_x_Excur','swing_x_Excur','stance_y_Excur','swing_y_Excur','stanceSlipFactor'};

%Put together the numeric data with text data containing headers, filename
ls_compiledCellArray = cell(1,12); %note this has one additional column for the filename information

%Create temporary cell for to hold column with fileName
c_ls_saveCell1 = cell(size(legs_strides_SaveMatrix,1), 1);
[c_ls_saveCell1{:,1}] = deal(filePrefix);

%Convert data from matrix to cell array for saving:
c_ls_saveCell2 =  num2cell(legs_strides_SaveMatrix);

%Compile filename column with data columns into a single cell array
c_saveCell = horzcat(c_ls_saveCell1,c_ls_saveCell2);

%Create a headers row for the top of the spreadsheet:
[ls_compiledCellArray{1,1:end}] = deal(legs_strides_Header{:});
%Put the headers row together with the data:
ls_compiledCellArray = vertcat(ls_compiledCellArray,c_saveCell);

%Save the cell array to a file:
cell2csv([kPathname,filePrefix, '_Data_Legs_Strides.csv' ], ls_compiledCellArray)

% Calculate relative leg phases

legPhaseDiffs = nan(size(rotatedKineData,1),7);
for i= 1:8
    legPhaseDiffs(:,i) = rad2deg(refPhase(:,refLeg) - refPhase(:,i));
end

%Create a big output matrix with body motion and relative phases with
%respect to reference leg (L2)

meanBodyPt(:,1) = mean(filtRotKineData(:,xBodyIndex),2);
meanBodyPt(:,2) = mean(filtRotKineData(:,yBodyIndex),2);

[~,bodyVelXY] = SplineInterp_wSPAPS(meanBodyPt, time, kin_sptol, 1);

%atan2d(Y,X)
bodyOrientation_deltaX = filtRotKineData(:,xBodyIndex(end)) - filtRotKineData(:,xBodyIndex(1));
bodyOrientation_deltaY = filtRotKineData(:,yBodyIndex(end)) - filtRotKineData(:,yBodyIndex(1));
bodyYawAngle = unwrap(atan2d(bodyOrientation_deltaY,bodyOrientation_deltaX));
bodyYawAngle = bodyYawAngle - mean(bodyYawAngle);
bodyVelMag = sqrt(bodyVelXY(:,1).^2 + bodyVelXY(:,2).^2);
bodyVelAng = unwrap(atan2d(bodyVelXY(:,2),bodyVelXY(:,1)));

fbd = figure; 
hold on;
subplot(3,1,1)
plot(time,bodyYawAngle);
    ylabel('Body Yaw Angle (deg)')
subplot(3,1,2)
plot(time,bodyVelMag);
ylabel('Body Velocity (mm/s)')
subplot(3,1,3)
plot(time,bodyVelAng);
ylabel('Body Velocity Angle (deg)')
xlabel('Time (s)')
subplot(3,1,1)
 hold on;
title([filePrefix ': Body Trajectory Diagram'] )
[~] = SaveFigAsPDF(fbd,kPathname,filePrefix,'_BodyTrajectory');

RefLegFootOffIndex =  combinedEvents{refLeg}(combinedEvents{refLeg}(:,2)==-1,1);

strideLength = nan(length(RefLegFootOffIndex),1);
strideVelocity = nan(length(RefLegFootOffIndex),1);
strideLength(2:end) = sqrt((diff(meanBodyPt(RefLegFootOffIndex,1)).^2 + diff(meanBodyPt(RefLegFootOffIndex,2)).^2));%in mm
strideVelocity(2:end) = strideLength(2:end)./(diff(RefLegFootOffIndex)./framerate); %mm per second

RefLegFootOffTimes = RefLegFootOffIndex./framerate;

body_phase_SaveMatrix = [[1:length(RefLegFootOffTimes)]' RefLegFootOffTimes bodyVelMag(RefLegFootOffIndex) bodyVelAng(RefLegFootOffIndex) bodyYawAngle(RefLegFootOffIndex) strideLength strideVelocity legPhaseDiffs(RefLegFootOffIndex,:)];
body_phase_Headers = {'fileName', 'strideNumber','RefLegFootOffTimes','bodyVelMag','bodyVelAng','bodyYawAngle','strideLength','strideVelocity','legPhaseDiff1','legPhaseDiff2', 'legPhaseDiff3','legPhaseDiff4','legPhaseDiff5','legPhaseDiff6','legPhaseDiff7','legPhaseDiff8'};

%Put together the numeric data with text data containing headers, filename
bp_compiledCellArray = cell(1,16);
%Create temporary cell for to hold column with fileName
c_bp_saveCell1 = cell(size(body_phase_SaveMatrix,1), 1);
[c_bp_saveCell1{:,1}] = deal(filePrefix);
%Convert data from matrix to cell array for saving:
c_bp_saveCell2 =  num2cell(body_phase_SaveMatrix);
%Compile filename column with data columns into a single cell array
c_bp_saveCell = horzcat(c_bp_saveCell1,c_bp_saveCell2);
%Create a headers row for the top of the spreadsheet:
[bp_compiledCellArray{1,1:end}] = deal(body_phase_Headers{:});
%Put the headers row together with the data:
bp_compiledCellArray = vertcat(bp_compiledCellArray,c_bp_saveCell);

%Save the cell array to a file:
cell2csv([kPathname,filePrefix, '_Data_Body_Phase.csv' ], bp_compiledCellArray)

out_LS_compiledArray = ls_compiledCellArray(2:end,:);
out_bp_compiledArray = bp_compiledCellArray(2:end,:);

% % plot & save bodyVelMag with locations of RefLegFootOffIndex marked, for possible
% % debugging/interest later
% f14 = figure;
% plot(time,bodyVelMag);
% hold on
% plot (time(RefLegFootOffIndex),bodyVelMag(RefLegFootOffIndex),'X','MarkerSize',12);
% xlabel('Time (s)')
% ylabel('COM velocity (mm/s)')
% legend ('Velocity','RefLegFootOffIndex')
% title([filePrefix ': COM Velocity'] )
% [~] = SaveFigAsPDF(f14,kPathname,filePrefix,'_BodyVelMag_RefLegFOInd');
% close(f14);

end


function [legData, xNewLegs, yNewLegs] = PullOutLegCoords (kineData,xLegsIndex,yLegsIndex)
% PullOutLegCoords takes leg coordinates from imported CSV data from
% ProAnalyst, and puts it in a matrix.

% create empty dataset & new XYs (what data type is this??)
legData = nan(size(kineData,1),length(xLegsIndex).*2);
xNewLegs = [1:2:size(legData,2)];
yNewLegs =  [2:2:size(legData,2)];

% pull out leg XY data
% in data, all rows with column numbers in xNew - fill with stuff in
% tempData column numbers in xPts
legData(:,xNewLegs) = kineData(:,xLegsIndex);
legData(:,yNewLegs) = kineData(:,yLegsIndex);
end

function [figH] = SaveFigAsPDF(figH,defDir,baseFNameString,suffixString)

if ~exist('suffixString','var')
    suffixString = [];
end
set(figH,'units', 'normalized'); set(figH,'Position', [0 0.0364583 1 0.875]);
figFilename = [defDir baseFNameString suffixString '.pdf'];
saveas(figH,figFilename,'pdf');
end


function [legLengths,leg_Y_RelToCom] = CalcLegLengths (legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols)

legLengths = nan(nRows,length(filtXLegCols));
leg_Y_RelToCom = nan(nRows,length(filtXLegCols));

for i=1:8
    % calculate leg lengths
    legLengths(:,i) = sqrt((legDataRelToCentroid(:,filtXLegCols(i)).^2)+(legDataRelToCentroid(:,filtYLegCols(i)).^2));
    %subtract mean leg lengths
    leg_Y_RelToCom(:,i) = legLengths(:,i) - nanmean(legLengths(:,i));
end

end


function [legAngles,leg_X_RelToCom] = CalcLegAngles (legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols)

legAngles = nan(nRows, length(filtXLegCols));
leg_X_RelToCom = nan(nRows,length(filtXLegCols));

for i=1:8
    % calculate leg angles
    legAngles(:,i) = atan2d(legDataRelToCentroid(:,filtYLegCols(i)),legDataRelToCentroid(:,filtXLegCols(i)));
end
legAnglesRad = (pi/180).*(legAngles);
legAngles = rad2deg(unwrap(legAnglesRad));

for i=1:8
    %subtract mean leg angle
    leg_X_RelToCom(:,i) = legAngles(:,i) - nanmean(legAngles(:,i));
end


end

function [Rad] = deg2rad(Deg) %#ok<DEFNU>
% Covert an angle in radians to degrees

Rad = (pi/180) * Deg;

end


function [combinedEvents]= CombineHilbertVelocityEvents(hilbertSorted,vthreshSorted,c_footVel)

%Delete any hilbert events that are less than some threshold apart
h_diff = [NaN; diff(hilbertSorted)];
%Note the value here is arbitrary and may need to be adjusted:
%Currently set to 5 frames (5ms for spider data sampled at 1000Hz)
hilbertSorted(h_diff<5)=[];

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

