function [out_bp_compiledArray,body_phase_Headers, out_LS_compiledArray, legs_strides_Header, filePrefix, combinedEvents] = ProcessSpiderData(fileName,pathName)
% A script to process spider gait data (initially orb-weaver, intact trials)

%Utility variable
if isempty(strfind(computer,'MAC'))
    sep = '\';
else sep = '/';
end

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
    mkdir([kPathname '_Results' sep ]);
else
    kFilename = fileName;
    kPathname = pathName;
end
filePrefix = kFilename(1:(end-13));

%Create a new folder to save the results files into
savePathN = dirCheck([kPathname '_Results' sep ]);

%Check for previously processed file, to skip smoothing steps if possible
if ~exist([kPathname filePrefix '.mat'],'file')
    
    % set framerate
    framerate = 1000;
    %Number of legs in analysis
    numLegs = 8;    
    % Reference leg = L2 (Cyclical leg but not adjacent to a later ablated leg
    % or part of the same segment pair)
    refLeg = 2;
    %Leg ID labels for graphs:
    % In order digitised:
    leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};
    %For graphs in anatomical arrangment
    leg_labels_anat = {'R4','L4','R3','L3','R2','L2','R1','L1'};

    % load in data file (tempData) - dlmread / csvread
    tempData = csvread([kPathname kFilename]);
             
    %Checks for NaN values at the beginning and end
    % of the file and delete the corresponding rows.
    [tempData] = RemoveNaNRowsAtStartAndEnd(tempData);
        
    % Pull out time/frame vectors
    time = tempData(:,2);
    frames = tempData(:,1); %#ok<NASGU>
    
    % Get number of rows
    nRows = size(tempData,1);
    
    % Pull out kinematic data (separatefrom time/frame columns
    kineCoords = [3:size(tempData,2)];
    kineData = nan(nRows,length(kineCoords));
    newKineCoords = [1:size(kineData,2)];
    kineData(:,newKineCoords) = tempData(:,kineCoords);
    
    % Rotate the data to be in the body frame of reference, with forward
    % motion in the positive X direction       
    [rotatedKineData] = RotateDataToBodyFrame(kineData,[1]);
    
    vel_sptol = 0.001; % sptol for smoothing velocities      
    [tempSmoothData,kineVel] = SplineInterp_wSPAPS(rotatedKineData, time, vel_sptol, 1);
    x_Cols = [1:2:size(rotatedKineData,2)]; % all X coords
    y_Cols = [2:2:size(rotatedKineData,2)]; % all Y coords

    % plot some legs raw resultant vel and smoothed resultant vel to check sptol
    f5=figure;
    plot(rotatedKineData(:,x_Cols),rotatedKineData(:,y_Cols),'r');
    hold on;
    plot(tempSmoothData(:,x_Cols), tempSmoothData(:,y_Cols),'b');
    xlabel('x-coordinates')
    ylabel('y-coordinates')
    title([filePrefix ': Kinematic data: red = raw, blue = filtered | sptol = ' num2str(vel_sptol)])
    [~] = SaveFigAsPDF(f5,savePathN,filePrefix,'_FiltKineData');
    close(f5);
    
    rotatedKineData = tempSmoothData;
    clear tempData tempSmoothData;
 
    %Get leg and body column IDs depending on file type case
    %Also adds NaN (empty) 'placeholder' columns for the 
    %Trials with leg ablations, so the gait diagrams are consistently
    %Based on the 'original' leg configuration
    [~,~,~,kineData] = GetColumnIndicesBasedOnFileName(kineData,kFilename);
    [~,~,~,kineVel] = GetColumnIndicesBasedOnFileName(kineVel,kFilename);
    [bodyPts,legPts,rLpts,rotatedKineData,conditionStr] = GetColumnIndicesBasedOnFileName(rotatedKineData,kFilename);
        
    x_Cols = [1:2:size(rotatedKineData,2)]; % all X coords
    y_Cols = [2:2:size(rotatedKineData,2)]; % all Y coords
    
         
    %Plot the original and rotated data for 'reality check'
    f1=figure;
    plot(kineData(:,x_Cols), kineData(:,y_Cols),'r');
    hold on;
    plot(kineData(1,x_Cols), kineData(1,y_Cols),'ro');
    plot(rotatedKineData(:,x_Cols), rotatedKineData(:,y_Cols),'b');
    plot(rotatedKineData(1,x_Cols), rotatedKineData(1,y_Cols),'bo');
    xlabel('x-coordinates')
    ylabel('y-coordinates')
    title([filePrefix ': Rotated kinematic data: red= raw, blue, rotated'])
    [~] = SaveFigAsPDF(f1,savePathN,filePrefix,'_RotatedData');
    close(f1);
    
    % Leg velocity data:
     smLegXVel = nan(size(rotatedKineData,1),length(legPts));
     smLegYVel = nan(size(rotatedKineData,1),length(legPts));
       
    smLegXVel(:,rLpts)= kineVel(:,x_Cols(legPts(rLpts)));    
    smLegYVel(:,rLpts)= kineVel(:,y_Cols(legPts(rLpts)));
    % calculate  resultant velocity (smTotVel) from smoothed XY velocities (smVel)
        smLegTotVel = sqrt(smLegXVel.^2 + smLegYVel.^2);
               
    %Calculate some velocity values for threshold detection of foot contact
    mean_legVel = nanmean(abs(smLegTotVel));
    std_legVel = nanstd(abs(smLegTotVel));
    t_velThreshNum = 0.5; % default velocity threshold number
    
        legContacts = zeros(size(smLegTotVel,1),length(leg_labels));
        legVel_InContact = nan(size(smLegTotVel,1),length(leg_labels));
        velThreshold = mean_legVel - t_velThreshNum.*(std_legVel);
        
        % For each foot, detect foot contact periods
        for i = 1:numLegs
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
        [~] = SaveFigAsPDF(f3,savePathN,filePrefix,'_VelocityDetection');
        
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
        ylim([1,numLegs])
        title( [filePrefix ': Initial Gait diagram: Foot velocity only.'])
        set(gca,'YTick',[1 2 3 4 5 6 7 8])
        set(gca,'YTickLabel',leg_labels_anat)
        [~] = SaveFigAsPDF(f4,savePathN,filePrefix,'_GaitDiagram1');
         
    % pull out just leg data for calculating lengths and angles
    [filtRLegData, filtXLegCols, filtYLegCols] = PullOutLegCoords(rotatedKineData,x_Cols,y_Cols,legPts);
    %[legData, xNewLegs, yNewLegs] = PullOutLegCoords(kineData,x_Cols,y_Cols,legPts)
    
    % Calculate legs Relative to body (Leg X/Y - COM X/Y)
    legDataRelToCOM = nan(nRows,size(filtRLegData,2));
    for i=1:numLegs
        legDataRelToCOM(:,filtXLegCols(i)) = filtRLegData(:,filtXLegCols(i)) - nanmean(rotatedKineData(:,x_Cols(bodyPts)),2);
        legDataRelToCOM(:,filtYLegCols(i)) = filtRLegData(:,filtYLegCols(i)) - nanmean(rotatedKineData(:,y_Cols(bodyPts)),2);
    end
    
    % % calculate leg lengths & angles
    % [legLengths] = CalcLegLengths(legDataRelToCOM,nRows,filtXLegCols,filtYLegCols);
    % [legAngles] = CalcLegAngles(legDataRelToCOM,nRows,filtXLegCols,filtYLegCols);
    
    aveLegX = nan(1,length(filtXLegCols));
    aveLegY = nan(1,length(filtYLegCols));
    aveLegLength = nan(1,length(filtYLegCols));
    
    % Calculate legs Relative to body (Leg X/Y - COM X/Y)
    legDataRelToCentroid = nan(nRows,size(filtRLegData,2));
    for i=1:numLegs
        aveLegX(i) = mean(legDataRelToCOM(:,filtXLegCols(i)));
        aveLegY(i) = mean(legDataRelToCOM(:,filtYLegCols(i)));
        
        legDataRelToCentroid(:,filtXLegCols(i)) = legDataRelToCOM(:,filtXLegCols(i)) - aveLegX(i);
        legDataRelToCentroid(:,filtYLegCols(i)) = legDataRelToCOM(:,filtYLegCols(i)) - aveLegY(i);
        aveLegLength(i) = sqrt((mean(legDataRelToCOM(:,filtXLegCols(i)))).^2 +(mean(legDataRelToCOM(:,filtYLegCols(i)))).^2);
    end
    
    leg_X_RelToCom = legDataRelToCentroid(:,filtXLegCols);
    leg_Y_RelToCom = legDataRelToCentroid(:,filtYLegCols);
    
    legDataRelToCentroid(:,filtXLegCols) = legDataRelToCentroid(:,filtXLegCols)+ aveLegX(refLeg);
    legDataRelToCentroid(:,filtYLegCols) = legDataRelToCentroid(:,filtYLegCols)+ aveLegY(refLeg);
    
    % calculate leg lengths & angles
    [legLengthsCentroid,legLengthsMeanSub] = CalcLegLengths(legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols,numLegs); %#ok<ASGLU>
    [legAnglesCentroid,legAnglesMeanSub] = CalcLegAngles(legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols,numLegs); %#ok<ASGLU>
    
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
[~] = SaveFigAsPDF(flx,savePathN,filePrefix,'_Foot_XDisp');
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
[~] = SaveFigAsPDF(flx,savePathN,filePrefix,'_Foot_YDisp');
%close(flx);

%Plot 'overhead' view plot in cartesian coordinates
% Replaces previous polar plot
f10=(figure); %#ok<NASGU>
plot(legDataRelToCOM(:,filtXLegCols),legDataRelToCOM(:,filtYLegCols));
legend(leg_labels,'Location', 'eastoutside');
title([filePrefix ': Foot trajectory (rel to COM)']);
ylabel('y-coordinates');
xlabel('x-coordinates');
print([savePathN,filePrefix,'_FootTrajectoryPlot'],'-dpdf');
%close(f10);

meanBodyPt(:,1) = mean(rotatedKineData(:,x_Cols(bodyPts)),2);
meanBodyPt(:,2) = mean(rotatedKineData(:,y_Cols(bodyPts)),2);

[~,bodyVelXY] = SplineInterp_wSPAPS(meanBodyPt, time, vel_sptol, 1);
bodyVelMag = sqrt(bodyVelXY(:,1).^2 + bodyVelXY(:,2).^2);
bodyVelAng = unwrap(atan2d(bodyVelXY(:,2),bodyVelXY(:,1)));

%atan2d(Y,X)
bodyOrientation_deltaX = rotatedKineData(:,x_Cols(bodyPts(end))) - rotatedKineData(:,x_Cols(bodyPts(1)));
bodyOrientation_deltaY = rotatedKineData(:,y_Cols(bodyPts(end))) - rotatedKineData(:,y_Cols(bodyPts(1)));

bodyYawAngle = unwrap(atan2d(bodyOrientation_deltaY,bodyOrientation_deltaX));
bodyYawAngle = bodyYawAngle - mean(bodyYawAngle);

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
[~] = SaveFigAsPDF(fbd,savePathN,filePrefix,'_BodyTrajectory');

% leg_disp_RelToCoM = sqrt(leg_X_RelToCom.^2 + leg_Y_RelToCom.^2);
% legAngleRelToCoM = atan2d(leg_Y_RelToCom,leg_X_RelToCom);
% legBehindCoM = abs(legAngleRelToCoM)>90;
% leg_disp_RelToCoM(legBehindCoM) = leg_disp_RelToCoM(legBehindCoM).*-1;
% 
% figure; 
% plot(time,leg_disp_RelToCoM(:,2),'k')
% hold on;
% plot(time,leg_X_RelToCom(:,2),'r')
% plot(time,leg_Y_RelToCom(:,2),'b')

% calculate Hilbert phases (normal and inverted)
for i=1:numLegs
    hilbert_phase(:,i) = angle(hilbert(leg_X_RelToCom(:,i)));
    hilbert_phase_inverted(:,i) = angle(hilbert(leg_X_RelToCom(:,i)).*-1);
%     hilbert_phase(:,i) = angle(hilbert(legAnglesMeanSub(:,i)));
%     hilbert_phase_inverted(:,i) = angle(hilbert(legAnglesMeanSub(:,i)).*-1);   
end

%These variables  need to be pre-allocated to ensure
% that they are specified to have the correct size/dimensions
footOff_H1 = cell(numLegs,1);
footDown_H2 = cell(numLegs,1);
footDown_I = cell(numLegs,1);
footOff_I = cell(numLegs,1);

%Add body velocity, velocity angle and body yaw to these output variables
% To calculate for every leg. 
stridePeriod = cell(numLegs,1);
stancePeriod = cell(numLegs,1);
swingPeriod = cell(numLegs,1);
dutyFactor = cell(numLegs,1);

strideAveVel = cell(numLegs,1);
strideAveVelAng = cell(numLegs,1);
strideAveYawAng = cell(numLegs,1);
strideDeltaVel = cell(numLegs,1);
strideDeltaVelAng = cell(numLegs,1);
strideDeltaYawAng = cell(numLegs,1);

stance_x_Excur = cell(numLegs,1);
swing_x_Excur = cell(numLegs,1);
stance_y_Excur = cell(numLegs,1);
swing_y_Excur = cell(numLegs,1);
stanceSlipFactor = cell(numLegs,1);

refPhase = NaN(size(hilbert_phase,1),numLegs);

% set up for hilbert event detection while loop for user input
hilEventDet_thresh = 0.06; % default threshold for detecting events from Hilbert

    %Need to recreate these variables in each iteration of while loop
    %To ensure they are recalculated correctly if the hilbert threshold is adjusted
    newGaitDiagram = nan(size(legContacts));
    hilbertEvents = cell(numLegs,1);
    velocityEvents = cell(numLegs,1);
    combinedEvents = cell(numLegs,1);
    
    for i=1:numLegs
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
            
            %Update gait diagram
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
            
            % After updating raw gait diagram, delete any fD events that
            % occur before foot-off events, because we are calculating
            % Values only for 'full stride cycles', defined between
            % foot-off and the next foot-off
            if fD_Events(1)< fO_Events(1)           
                 fD_Events(1) = [];
                 delFD = find(newEventsSorted(:,2)== 1,1,'first');
                 newEventsSorted(delFD,:)=[];                
            end
                
           if fD_Events(end)>fO_Events(end)
              fD_Events(end) = [];
              delFD = find(newEventsSorted(:,2)== 1,1,'last');
              newEventsSorted(delFD,:)=[];                
           end

            % Do stride based calculations for all full stride cycles 
            % Strides are calculated from foot off to the next foot off
            % Index to adjacent '-1' events, which correspond to foot off
            stridePeriod{i} = [diff(fO_Events)]./framerate;
            numStrides = length(stridePeriod{i,1});

            %Save event information, save to CSV file later
            combinedEvents{i} = newEventsSorted;
                            
            halfPeriods_combined = diff(newEventsSorted,1);           
            stancePeriod{i} = (halfPeriods_combined(halfPeriods_combined(:,2)==-2,1))./framerate;
            swingPeriod{i} = (halfPeriods_combined(halfPeriods_combined(:,2)==2,1))./framerate;
            t_end = min([length(stridePeriod{i}) length(stancePeriod{i})]);
            dutyFactor{i} =  stancePeriod{i}(1:t_end)./ stridePeriod{i}(1:t_end);
                        
            for S = 1:numStrides              
                    strideAveVel{i,1}(S) = mean(bodyVelMag(fO_Events(S,1):fO_Events(S+1,1),1));
                    strideAveVelAng{i,1}(S) = mean(bodyVelAng(fO_Events(S,1):fO_Events(S+1,1),1));
                    strideAveYawAng{i,1}(S) = mean(bodyYawAngle(fO_Events(S,1):fO_Events(S+1,1),1));
                    strideDeltaVel{i,1}(S) = max(bodyVelMag(fO_Events(S,1):fO_Events(S+1,1),1)) - min(bodyVelMag(fO_Events(S,1):fO_Events(S+1,1),1));
                    strideDeltaVelAng{i,1}(S) = max(bodyVelAng(fO_Events(S,1):fO_Events(S+1,1),1)) - min(bodyVelAng(fO_Events(S,1):fO_Events(S+1,1),1));
                    strideDeltaYawAng{i,1}(S) = max(bodyYawAngle(fO_Events(S,1):fO_Events(S+1,1),1)) - min(bodyYawAngle(fO_Events(S,1):fO_Events(S+1,1),1));
            end
            
            for S = 1:numStrides
                %Calculate stance and swing averages within each stride
                %Find the current foot down event
                %Calculations only done for complete stride cycles
                %So some stance and swing phases will be lost, 
                % if they occur outside the complete stride cycles
                % (as defined above, from foot off to the next foot off
                %That's OK, for consistency with the calculations
                c_FDevent = fD_Events(fD_Events(:,1)>fO_Events(S,1) & fD_Events(:,1)<fO_Events(S+1,1),1);                
            
                    c_foot_X_stance = leg_X_RelToCom(c_FDevent:fO_Events(S+1,1),i);
                    c_foot_Y_stance = leg_Y_RelToCom(c_FDevent:fO_Events(S+1,1),i);
                    %Calculate 'slip factor'
                    c_stanceSlipFactor = sum(~isnan(gaitDiagramData(c_FDevent:fO_Events(S+1,1),i)))./sum(~isnan(newGaitDiagram(c_FDevent:fO_Events(S+1,1),i)));
                    c_foot_X_swing = leg_X_RelToCom(fO_Events(S,1):c_FDevent,i);
                    c_foot_Y_swing = leg_Y_RelToCom(fO_Events(S,1):c_FDevent,i);
                
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
                       
            % find values where dutyFactor is <0 - this indicates some
            % problematic stride/stance periods. Replace these strides with NaN
            badValues = find([dutyFactor{i,:}]<0 | [dutyFactor{i,:}]>1);
            % (Note that the index badValues comes after the cell index.)
            dutyFactor{i,1}(badValues) = NaN;
            stancePeriod{i,1}(badValues) = NaN;
            swingPeriod{i,1}(badValues) = NaN;
            stridePeriod{i,1}(badValues) = NaN;
            
            strideAveVel{i,1}(badValues) =  NaN;
            strideAveVelAng{i,1}(badValues) =  NaN;
            strideAveYawAng{i,1}(badValues) =  NaN;
            strideDeltaVel{i,1}(badValues) =  NaN;
            strideDeltaVelAng{i,1}(badValues) =  NaN;
            strideDeltaYawAng{i,1}(badValues) =  NaN;

            swing_x_Excur{i,1}(badValues) =  NaN;
            swing_y_Excur{i,1}(badValues)=  NaN;
            stance_x_Excur{i,1}(badValues) =  NaN;
            stance_y_Excur{i,1}(badValues) =  NaN;
            stanceSlipFactor{i,1}(badValues) =  NaN;
            
            numStrides = length(stridePeriod{i,1}); %#ok<NASGU>

        else %If there are not enough events to calculate full stride cycles for this leg,  enter NaNs
            newEventsSorted = NaN;
            combinedEvents{i} = newEventsSorted;
            stancePeriod{i} = NaN;
            swingPeriod{i} = NaN;
            stridePeriod{i} = NaN;
            dutyFactor{i} =  NaN ;
            
            strideAveVel{i} =  NaN;
            strideAveVelAng{i} =  NaN;
            strideAveYawAng{i} =  NaN;
            strideDeltaVel{i} =  NaN;
            strideDeltaVelAng{i} =  NaN;
            strideDeltaYawAng{i} =  NaN;

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
    title([filePrefix ': Reference phases (with foot off/down events) | HilbertEvent threshold = ' num2str(hilEventDet_thresh)])
    [~] = SaveFigAsPDF(f14,savePathN,filePrefix,'_refPhase');
    close(f14);
  
    
% %plot reference phases
%     f14 = figure;
%     for i = 1:4
%         subplot(4,1,i)
%         hold on;
%         col_i = SpidColors(i,:);
%         col_i2 = SpidColors(i+4,:);
%         plot(time,hilbert_phase(:,i),'Color',col_i);
%         plot(time,hilbert_phase(:,i+4),'Color',col_i2);
%         % plot events for left legs (H2 = foot off)
%         plot(time(footDown_H2{i}), hilbert_phase(footDown_H2{i},i),'ko'); % foot off
%         plot(time(footOff_H1{i}), hilbert_phase(footOff_H1{i},i),'kX'); %foot on
%         % plot events for right legs (H1 = foot off)
%         plot(time(footOff_H1{i+4}), hilbert_phase(footOff_H1{i+4},i+4),'ko'); % foot off
%         plot(time(footDown_H2{i+4}), hilbert_phase(footDown_H2{i+4},i+4),'kX'); % foot on
%         
%         legend(['L' num2str(i)],['R' num2str(i)],'FO', 'FD');
%         if i == 2
%             ylabel('Hilbert phase')
%         end
%     end
%     xlabel('Time (s)')
%     subplot(4,1,1)
%     hold on;
%     title([filePrefix ': Hilbert phases (with foot off/down events) | HilbertEvent threshold = ' num2str(hilEventDet_thresh)])
%     [~] = SaveFigAsPDF(f14,savePathN,filePrefix,'_HilbertPhase');
%     %close(f14);
    
%     %plot reference phases
%     f14 = figure;
%     for i = 1:4
%         subplot(4,1,i)
%         hold on;
%         col_i = SpidColors(i,:);
%         col_i2 = SpidColors(i+4,:);
%         plot(time,hilbert_phase_inverted(:,i),'Color',col_i);
%         plot(time,hilbert_phase_inverted(:,i+4),'Color',col_i2);
%         % plot events for left legs (H2 = foot off)
%         plot(time(footDown_H2{i}), hilbert_phase_inverted(footDown_H2{i},i),'ko'); % foot off
%         plot(time(footOff_H1{i}), hilbert_phase_inverted(footOff_H1{i},i),'kX'); %foot on
%         % plot events for right legs (H1 = foot off)
%         plot(time(footOff_H1{i+4}), hilbert_phase_inverted(footOff_H1{i+4},i+4),'ko'); % foot off
%         plot(time(footDown_H2{i+4}), hilbert_phase_inverted(footDown_H2{i+4},i+4),'kX'); % foot on
%         
%         legend(['L' num2str(i)],['R' num2str(i)],'FO', 'FD');
%         if i == 2
%             ylabel('Hilbert phase')
%         end
%     end
%     xlabel('Time (s)')
%     subplot(4,1,1)
%     hold on;
%     title([filePrefix ': Hilbert Inverted phases (with foot off/down events) | HilbertEvent threshold = ' num2str(hilEventDet_thresh)])
%     [~] = SaveFigAsPDF(f14,savePathN,filePrefix,'_HilbertInverted');
%    % close(f14);    
    
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
    [~] = SaveFigAsPDF(f12,savePathN,filePrefix,'_GaitDiagram_Combined');
    
% SAVE OUT STRIDE PARAMETERS AS CSV - LEGS/STRIDES

%Create empty matrix for final leg matrix (all legs compiled)
legs_strides_SaveMatrix = [];

for i=1:numLegs %For all legs
    
    cSeg = ceil(i/2)-1;    
    %Right = 0,  Left = 1
    r_v_L = ceil(mod(i,2)/2);
    
    c_nRws = max([length(stancePeriod{i}) length(swingPeriod{i}) length(stridePeriod{i})]);
    %Matrix for current leg data
    legs_strides_tempMatrix = nan(c_nRws,19); %create temp matrix full of NaNs
    
    legs_strides_tempMatrix(1:c_nRws,1) = repmat(i,c_nRws,1); %leg number
    %    legs_strides_tempMatrix(1:c_nRws,1) = repmat(i,c_nRws,1); %leg number

    legs_strides_tempMatrix(1:c_nRws,2) = [1:c_nRws]; %stride number
    
    legs_strides_tempMatrix(1:(length(strideAveVel{i})),3) = strideAveVel{i}; %note indexing here in case this variable has less than the max strides
    legs_strides_tempMatrix(1:(length(strideAveVelAng{i})),4) = strideAveVelAng{i}; 
    legs_strides_tempMatrix(1:(length(strideAveYawAng{i})),5) = strideAveYawAng{i}; 
    legs_strides_tempMatrix(1:(length(strideDeltaVel{i})),6) = strideDeltaVel{i}; 
    legs_strides_tempMatrix(1:(length(strideDeltaVelAng{i})),7) = strideDeltaVelAng{i}; 
    legs_strides_tempMatrix(1:(length(strideDeltaYawAng{i})),8) = strideDeltaYawAng{i}; 
    legs_strides_tempMatrix(1:(length(stridePeriod{i})),9) = stridePeriod{i}; 
    legs_strides_tempMatrix(1:(length(swingPeriod{i})),10) = swingPeriod{i};
    legs_strides_tempMatrix(1:(length(stancePeriod{i})),11) = stancePeriod{i};
    legs_strides_tempMatrix(1:(length(dutyFactor{i})),12) = dutyFactor{i};
    legs_strides_tempMatrix(1:(length(stance_x_Excur{i})),13) = stance_x_Excur{i};
    legs_strides_tempMatrix(1:(length(swing_x_Excur{i})),14) = swing_x_Excur{i};
    legs_strides_tempMatrix(1:(length(stance_y_Excur{i})),15) = stance_y_Excur{i};
    legs_strides_tempMatrix(1:(length(swing_y_Excur{i})),16) = swing_y_Excur{i};
    legs_strides_tempMatrix(1:(length(stanceSlipFactor{i})),17) = stanceSlipFactor{i};
    legs_strides_tempMatrix(1:c_nRws,18) = repmat(cSeg,c_nRws,1); 
    legs_strides_tempMatrix(1:c_nRws,19) = repmat(r_v_L,c_nRws,1); 

    %concatenate current leg data into the final leg matrix
    legs_strides_SaveMatrix = [legs_strides_SaveMatrix; legs_strides_tempMatrix];
end

%Create a cell array to save as a CSV that includes filename information and headers within the file.
%Filename is repeated in each row so that multiple files can be compiled into a single spreadsheet later.

%Create output header row
legs_strides_Header = {'fileName','Species','SubjectNum', 'AblationCond','Terrain','EggsacCond', 'legNumber', 'strideNumber', 'strideAveVel','strideAveVelAng','strideAveYawAng',  ...
    'strideDeltaVel','strideDeltaVelAng','strideDeltaYawAng', 'stridePeriod', 'swingPeriod', ...
    'stancePeriod','dutyFactor','stance_x_Excur','swing_x_Excur','stance_y_Excur','swing_y_Excur','stanceSlipFactor','SegNum','R/L'};

%Put together the numeric data with text data containing headers, filename
ls_compiledCellArray = cell(1,length(legs_strides_tempMatrix)+6);
%note this has one additional column for the filename information

%Create temporary cell for to hold column with fileName and condition
%string
c_ls_saveCell1 = cell(size(legs_strides_SaveMatrix,1), 6);
[c_ls_saveCell1{:,1}] = deal(filePrefix);
[c_ls_saveCell1{:,2}] = deal(conditionStr{1}); %species
[c_ls_saveCell1{:,3}] = deal(conditionStr{2}); %SubjectNum
[c_ls_saveCell1{:,4}] = deal(conditionStr{3}); %AblationCond
[c_ls_saveCell1{:,5}] = deal(conditionStr{4}); %Terrain
[c_ls_saveCell1{:,6}] = deal(conditionStr{5}); %EggsacCond

%Convert data from matrix to cell array for saving:
c_ls_saveCell2 =  num2cell(legs_strides_SaveMatrix);

%Compile filename column with data columns into a single cell array
c_saveCell = horzcat(c_ls_saveCell1,c_ls_saveCell2);

%Create a headers row for the top of the spreadsheet:
[ls_compiledCellArray{1,1:end}] = deal(legs_strides_Header{:});
%Put the headers row together with the data:
ls_compiledCellArray = vertcat(ls_compiledCellArray,c_saveCell);

%Save the cell array to a file:
cell2csv([savePathN,filePrefix, '_Data_Legs_Strides.csv' ], ls_compiledCellArray);

% Calculate relative leg phases
unWr_refPhase = unwrap(refPhase);
%figure; plot(unWr_refPhase)
legPhaseDiffs = nan(size(rotatedKineData,1),7);
for i= 1:numLegs
    legPhaseDiffs(:,i) = rad2deg(unWr_refPhase(:,refLeg) - unWr_refPhase(:,i));   
end
%figure; plot(legPhaseDiffs)

     rL_nStrides = length(stridePeriod{refLeg,1});
     rl_Events = combinedEvents{refLeg};
      rf_fOEvents = rl_Events(rl_Events(:,2)== -1,1);
      
      aveLegPhaseDiff = nan(rL_nStrides+1,numLegs);
      for S = 1:rL_nStrides
          aveLegPhaseDiff(S+1,:)= median(legPhaseDiffs(rf_fOEvents(S,1):rf_fOEvents(S+1,1),:));         
      end      
      aveLegPhaseDiff = mod(aveLegPhaseDiff,360);
      
%Create a big output matrix with body motion and relative leg phases with
%respect to reference leg (L2)

RefLegFootOffIndex =  combinedEvents{refLeg}(combinedEvents{refLeg}(:,2)==-1,1);

strideLength = nan(length(RefLegFootOffIndex),1);
strideLength(2:end) = sqrt((diff(meanBodyPt(RefLegFootOffIndex,1)).^2 + diff(meanBodyPt(RefLegFootOffIndex,2)).^2));%in mm

% strideVelocity = nan(length(RefLegFootOffIndex),1);
% strideVelocity(2:end) = strideLength(2:end)./(diff(RefLegFootOffIndex)./framerate); %mm per second

RefLegFootOffTimes = RefLegFootOffIndex./framerate;

body_phase_SaveMatrix = [[NaN, 1:(length(RefLegFootOffTimes)-1)]' RefLegFootOffTimes ...
    [NaN strideAveVel{refLeg}]' [NaN strideDeltaVel{refLeg}]'...
    [NaN strideAveVelAng{refLeg}]' [NaN strideDeltaVelAng{refLeg}]'...
    [NaN strideAveYawAng{refLeg}]' [NaN strideDeltaYawAng{refLeg}]'...
    strideLength [NaN;dutyFactor{refLeg}] [NaN; stanceSlipFactor{refLeg}] aveLegPhaseDiff];

%eventBasedLegPhaseDiff = mod(legPhaseDiffs(RefLegFootOffIndex,:),360); 

body_phase_Headers = {'fileName','Species','SubjectNum', 'AblationCond','Terrain','EggsacCond', 'strideNumber','RefLegFootOffTimes',...
    'strideAveVel', 'strideDeltaVel',...
    'strideAveVelAng','strideDeltaVelAng',...
    'strideAveYawAng', 'strideDeltaYawAng', ...
    'strideLength','dutyFactor', 'stanceSlipFactor',...
    'legPhaseDiffL1','legPhaseDiffL2', 'legPhaseDiffL3','legPhaseDiffL4','legPhaseDiffR1','legPhaseDiffR2','legPhaseDiffR3','legPhaseDiffR4'};

%Put together the numeric data with text data containing headers, filename
bp_compiledCellArray = cell(1,length(body_phase_SaveMatrix)+6);
%Create temporary cell for to hold column with fileName
c_bp_saveCell1 = cell(size(body_phase_SaveMatrix,1), 6);
[c_bp_saveCell1{:,1}] = deal(filePrefix);
[c_bp_saveCell1{:,2}] = deal(conditionStr{1});%species
[c_bp_saveCell1{:,3}] = deal(conditionStr{2});%SubjectNum
[c_bp_saveCell1{:,4}] = deal(conditionStr{3});%AblationCond
[c_bp_saveCell1{:,5}] = deal(conditionStr{4});%Terrain
[c_bp_saveCell1{:,6}] = deal(conditionStr{5});%EddsacCond

%Convert data from matrix to cell array for saving:
c_bp_saveCell2 =  num2cell(body_phase_SaveMatrix);
%Compile filename column with data columns into a single cell array
c_bp_saveCell = horzcat(c_bp_saveCell1,c_bp_saveCell2);
%Create a headers row for the top of the spreadsheet:
[bp_compiledCellArray{1,1:end}] = deal(body_phase_Headers{:});
%Put the headers row together with the data:
bp_compiledCellArray = vertcat(bp_compiledCellArray,c_bp_saveCell);

%Save the cell array to a file:
cell2csv([savePathN,filePrefix, '_Data_Body_Phase.csv' ], bp_compiledCellArray);

%Export just the data portion of the array, as it will be compiled
%With the filename in the batch processing script
out_LS_compiledArray = ls_compiledCellArray(2:end,:);
out_bp_compiledArray = bp_compiledCellArray(2:end,:);

%Save Mat file
close all;
    save([kPathname filePrefix]) ;  
end


function[rotatedKineData] = RotateDataToBodyFrame(kineData,bodyPts)
    
nCols = size(kineData,2);
nPts = fix(nCols./2);
x_Cols = [1:2:size(kineData,2)]; % all X coords
y_Cols = [2:2:size(kineData,2)]; % all Y coords

    startIdx = find(sum(isnan(kineData),2) == 0, 1,'first' );
    endIdx = find(sum(isnan(kineData),2) == 0, 1,'last' );
    %Set origin arbirarily to first body point in data
    originPt = kineData(startIdx,x_Cols(bodyPts(1)):y_Cols(bodyPts(1)));
    kineData = kineData - repmat(originPt,size(kineData,1),nPts);
    
    % Calculate the net direction of motion and a rotation matrix
    deltaXtravelled = mean(kineData(endIdx,x_Cols(bodyPts))- kineData(startIdx,x_Cols(bodyPts)));
    deltaYtravelled = mean(kineData(endIdx,y_Cols(bodyPts))- kineData(startIdx,y_Cols(bodyPts)));
    rotation_angle = (atan2(deltaYtravelled,deltaXtravelled)).*-1;
    
    rotation_angle_deg = rad2deg(rotation_angle); %#ok<NASGU>
    rotation_matrix =  [  cos(rotation_angle), -sin(rotation_angle);
        sin(rotation_angle),  cos(rotation_angle)];
    
    rotatedKineData = nan(size(kineData,1),size(kineData,2));
    
    %Rotate the data so that the primary direction of motion is fore-aft
    for k=1:nPts
        c_point= [x_Cols(k)  y_Cols(k)];
        c_vector = kineData(:,c_point)';
        newVector = rotation_matrix*c_vector;
        newVector = newVector';
        rotatedKineData(:,c_point) = newVector;
    end
    
    if  deltaXtravelled <0
        rotatedKineData(:,yPts) = rotatedKineData(:,yPts).*-1;
    end
    
    %Read just origin vertical minimum
    yMin = min(min(rotatedKineData(:,y_Cols)));
    rotatedKineData(:,y_Cols) = rotatedKineData(:,y_Cols) - repmat(yMin,length(rotatedKineData(:,y_Cols)),nPts);

% %Now do a 2nd correction based on shorter term direction of motion
% Note initial try of this looked crazy-  will have to do a moving
% average of the data to rotate over time scales of approximately 1 stride
% not the instantaneous changes in direction

% [smKineData] = SplineInterp_wSPAPS(kineData, time, 0.001, 1);
%     startIdx = find(sum(isnan(smKineData),2) == 0, 1,'first' );   
%     originPt = smKineData(startIdx,x_Cols(bodyPts)(1):y_Cols(bodyPts)(1)); %XY coords of bodyCOM
%     smKineData = smKineData - repmat(originPt,size(smKineData,1),nPts);
%     
%     % Calculate the rotation matrix
%         newRotatedKineData = nan(size(smKineData,1),size(smKineData,2));
%         newRotatedKineData(1,:) = rotatedKineData(1,:);
%         
%     for r_i = 2:size(smKineData,1)
%         
%         c_deltaX = mean(smKineData(r_i,x_Cols(bodyPts))- kineData(r_i -1,x_Cols(bodyPts)));
%         c_deltaY = mean(smKineData(r_i,y_Cols(bodyPts))- kineData(r_i -1,y_Cols(bodyPts)));
%         c_rotation_angle = (atan2(c_deltaY,c_deltaX)).*-1;
%         
%     rotation_matrix =  [  cos(c_rotation_angle), -sin(c_rotation_angle);
%         sin(c_rotation_angle),  cos(c_rotation_angle)];
%     %Rotate the data so that the body's direction of motion the x-axis
%     for k=1:nPts
%         c_point= [x_Cols(k)  y_Cols(k)];
%         c_vector = kineData(r_i,c_point)';
%         newVector = rotation_matrix*c_vector;
%         newVector = newVector';
%         newRotatedKineData(r_i,c_point) = newVector;
%     end
%     end
%     %Re-adjust origin vertical minimum
%     yMin = min(min(newRotatedKineData(:,y_Cols)));
%     newRotatedKineData(:,y_Cols) = newRotatedKineData(:,y_Cols) - repmat(yMin,length(newRotatedKineData(:,y_Cols)),nPts);
          
end

function [legData, xNewLegs, yNewLegs] = PullOutLegCoords(kineData,x_Cols,y_Cols,legPts)
% PullOutLegCoords takes leg coordinates from imported CSV data from
% ProAnalyst, and puts it in a matrix.

% create empty dataset & new XYs (what data type is this??)
legData = nan(size(kineData,1),length(x_Cols(legPts)).*2);
xNewLegs = [1:2:size(legData,2)];
yNewLegs =  [2:2:size(legData,2)];

% pull out leg XY data
% in data, all rows with column numbers in xNew - fill with stuff in
% tempData column numbers in xPts
legData(:,xNewLegs) = kineData(:,x_Cols(legPts));
legData(:,yNewLegs) = kineData(:,y_Cols(legPts));
end

function [figH] = SaveFigAsPDF(figH,defDir,baseFNameString,suffixString)

if ~exist('suffixString','var')
    suffixString = [];
end
set(figH,'units', 'normalized'); set(figH,'Position', [0 0.0364583 1 0.875]);
figFilename = [defDir baseFNameString suffixString '.pdf'];
saveas(figH,figFilename,'pdf');
end


function [legLengths,leg_Y_RelToCom] = CalcLegLengths(legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols,numLegs)

legLengths = nan(nRows,length(filtXLegCols));
leg_Y_RelToCom = nan(nRows,length(filtXLegCols));

for i=1:numLegs
    % calculate leg lengths
    legLengths(:,i) = sqrt((legDataRelToCentroid(:,filtXLegCols(i)).^2)+(legDataRelToCentroid(:,filtYLegCols(i)).^2));
    %subtract mean leg lengths
    leg_Y_RelToCom(:,i) = legLengths(:,i) - nanmean(legLengths(:,i));
end

end


function [legAngles,leg_X_RelToCom] = CalcLegAngles(legDataRelToCentroid,nRows,filtXLegCols,filtYLegCols,numLegs)

legAngles = nan(nRows, length(filtXLegCols));
leg_X_RelToCom = nan(nRows,length(filtXLegCols));

for i=1:numLegs
    % calculate leg angles
    legAngles(:,i) = atan2d(legDataRelToCentroid(:,filtYLegCols(i)),legDataRelToCentroid(:,filtXLegCols(i)));
end
legAnglesRad = (pi/180).*(legAngles);
legAngles = rad2deg(unwrap(legAnglesRad));

for i=1:numLegs
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

function [bodyPts,legPts,rLpts,kineData,conditionStr] = GetColumnIndicesBasedOnFileName(kineData,kFilename)
%This function assigns column indices for body and legs based on the file type cases listed below
%Specify column indices for body and legs

% Species Individual Condition - Species| SubjectNum | Ablation | Terrain | Eggsac
conditionStr = {kFilename(4:7) kFilename(8:11) kFilename(13:18) kFilename(20) kFilename(22:24)};

if strcmp(kFilename(1:2),'03')==1 || strcmp(kFilename(1:2),'04')==1
    %Intact Wolf Spider trials with no egg sac,
    % 1-2 = bodyCOM X,Y
    % 3-4 = bodyBack X,Y
    % 5-6 = bodyFront X,Y
    % 7-8 = L1 X,Y
    % 9-10 = L2 X,Y
    % 11-12 = L3 X,Y
    % 13-14 = L4
    % 15-16 = R1
    % 17-18 = R2
    % 19-20 = R3
    % 21-22 = R4
    
    bodyPts = [1:3];
    legPts = [4:11];
    rLpts = [1:8];
    
elseif strcmp(kFilename(1:2),'01')==1 || strcmp(kFilename(1:2),'02')==1
    % Wolf spider trials with egg sac
    % 1-2 = bodyCOM X,Y
    % 3-4 = bodyBack X,Y
    % 5-6 = bodyFront X,Y
    % 7-8 = eggsac
    % 9-10 L1 X,Y
    % 11-12 = L2 X,Y
    % 13-14 = L3 X,Y
    % 15-16 = L4
    % 17-18 = R1
    % 19-20 = R2
    % 21-22 = R3
    % 23-24 = R4
    
    bodyPts = [1:4];
    legPts = [5:12];
    rLpts = [1:8];
    
elseif strcmp(kFilename(1:3),'12')==1 || strcmp(kFilename(1:2),'07')==1 || strcmp(kFilename(1:2),'08')==1
    % Wolf spider trials R4 ablated flat/rough no eggsac
    %  leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};
    % R4ablation trials - same leg columns, diff no of legs
    % Columns (for R4ablation data):
    % 1-2 = bodyCOM X,Y
    % 3-4 = bodyBack X,Y
    % 5-6 = bodyFront X,Y
    % 7-8 = L1 X,Y
    % 9-10 = L2 X,Y
    % 11-12 = L3 X,Y
    % 13-14 = L4
    % 15-16 = R1
    % 17-18 = R2
    % 19-20 = R3
    tempKineData = nan(size(kineData,1),size(kineData,2)+2);
    %Leave column of NaNs at the appropriate place in the dataset for
    %ablated legs
    tempKineData(:,1:size(kineData,2)) = kineData(:,1:size(kineData,2));
    kineData = tempKineData;
    
    bodyPts = [1:3];
    legPts = [4:11];
    rLpts = [1:7];
    
elseif strcmp(kFilename(1:3),'13')==1 || strcmp(kFilename(1:2),'09')==1 || strcmp(kFilename(1:2),'10')==1
    % Aran/Wolf L3 ablation R4 missing trials
    % L3ablation trials - same leg columns, diff no of legs
    %    leg_labels = {'L1','L2','L3','L4','R1','R2','R3','R4'};
    % Columns (for L3ablation R4missing data):
    % 1-2 = bodyCOM X,Y
    % 3-4 = bodyBack X,Y
    % 5-6 = bodyFront X,Y
    % 7-8 = L1 X,Y
    % 9-10 = L2 X,Y
    % 11-12 = L4
    % 13-14 = R1
    % 15-16 = R2
    % 17-18 = R3
    tempKineData = nan(size(kineData,1),size(kineData,2)+4);
    tempKineData(:,1:10) = kineData(:,1:10); %body pts, L1, L2
    tempKineData(:,13:size(tempKineData,2)-2) = kineData(:,11:size(kineData,2)); %L4, to end
    kineData = tempKineData;
    
    bodyPts = [1:3];
    legPts = [4:11];
    rLpts = [1 2 4 5 6 7];
    
else
    %Intact trials (Aran/Wolf) - '11' '05' '06'
    % 1-2 = bodyCOM X,Y
    % 3-4 = bodyBack X,Y
    % 5-6 = bodyFront X,Y
    % 7-8 = L1 X,Y
    % 9-10 = L2 X,Y
    % 11-12 = L3 X,Y
    % 13-14 = L4
    % 15-16 = R1
    % 17-18 = R2
    % 19-20 = R3
    % 21-22 = R4
    bodyPts = [1:3];
    legPts = [4:11];
    rLpts = [1:8];
    
    
end

end

function [tempData] = RemoveNaNRowsAtStartAndEnd(tempData)
    %This function checks for NaN values at the beginning and end
    % of the file and deletes the corresponding rows.
    % It avoids removing NaNs that appear in the middle of the file
    % because these can be interpolated in the filtering steps
    
   nCols = size(tempData,2);

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
end
