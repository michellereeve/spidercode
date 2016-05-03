% import data CSV from ProAnalyst as a numeric matrix
[kFilename, kPathname] = uigetfile(['.csv'], 'Choose coordinate data text (csv) file');


mkdir([kPathname '_CorrectedLensDist']);
saveDir = dirCheck([kPathname '_CorrectedLensDist']);

tempData = csvread([kPathname kFilename]);

%[rawData] = RemoveNaNRowsAtStartAndEnd(tempData);
[rawData] = RemoveNaNRowsAtEnd(tempData);

% Pull out time/frame vectors
time = rawData(:,2);
frames = rawData(:,1);

rawData = rawData(:,3:end); % remove frame and time columns

nCols = size(rawData,2);
x_Cols = [1:2:size(rawData,2)]; % all X coords
y_Cols = [2:2:size(rawData,2)]; % all Y coords

bodyFrontX = rawData(:,x_Cols(3));
bodyFrontY = rawData(:,y_Cols(3));

bodyBackX = rawData(:,x_Cols(2));
bodyBackY = rawData(:,y_Cols(2));

%To remove distortion and parallax effects
%Calculate body front-back length distance
% rescale all data at each timepoint so that this distance is constant 
%(equal to length at reference frame)
 
dX_bodyLength = bodyFrontX - bodyBackX;
dY_bodyLength = bodyFrontY - bodyBackY;
bodyLength = sqrt(dX_bodyLength.^2 + dY_bodyLength.^2);

% smooth body length so that the calculated scale factor isn't really jumpy - this is what caused the big spikes in data
sptol = 0.05; % sptol for smoothing body length (this value okayed by MD)
[bodyLength] = SplineInterp_wSPAPS(bodyLength, time, sptol, 1,time, 1); % will generate plot
    

% refFrame for INTACT Aran trials:
% 0022 = 1077
% 0023 = 912
% 0026 = 862
% 0027 = 555
% 0028 = 1029
% 0029 = 826
% 0030 = 530
% 0034 = 1339
% 0036 = 841

if strcmp(kFilename(10:11),'22')==1
    refFrame = 1077;
    
elseif strcmp(kFilename(10:11),'23')==1
    refFrame = 912;
    
elseif strcmp(kFilename(10:11),'26')==1
    refFrame = 862;

elseif strcmp(kFilename(10:11),'27')==1
    refFrame = 555;
    
elseif strcmp(kFilename(10:11),'28')==1
    refFrame = 1029;
    
elseif strcmp(kFilename(10:11),'29')==1
    refFrame = 826;
    
elseif strcmp(kFilename(10:11),'30')==1
    refFrame = 530;
    
elseif strcmp(kFilename(10:11),'34')==1
    refFrame = 1339;
    
elseif strcmp(kFilename(10:11),'36')==1
    refFrame = 841;
    
else 
    refFrame = 0;
   
end

scaleFactor = bodyLength./bodyLength(refFrame+1); % data starts at 0 frames so row 1 = frame 0, therefore row 1078= frame 1077 etc.
 
for i =1:nCols
rawData_N(:,i) = rawData(:,i)./scaleFactor;
bodyLengthN = bodyLength./scaleFactor;
end
 
figure; 
subplot(3,1,1)
hold on; 
plot(time, rawData_N(:,x_Cols));
%legend(headerNames)
ylabel('x-coordinates')
subplot(3,1,2)
hold on; 
plot(time, rawData_N(:,y_Cols));
ylabel('y-coordinates')
subplot(3,1,3)
plot(time,bodyLengthN)
ylabel('body length (mm)')
xlabel('time (s)')

newData = horzcat(frames,time,rawData_N);
newData = num2cell(newData);

savePathN = dirCheck([kPathname]);
%Save the cell array to a file:
cell2csv([saveDir, kFilename(1:8),'1',kFilename(10:end)], newData); % replaces second '0' with '1' to indicate it's corrected

