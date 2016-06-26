function [slowData, medData, fastData] = SetSpeedBins (data, VEL)
% split speed into bins

% 'MEDIUM' Speed: median ± 1IQR
medmin = median(VEL)-(0.5*iqr(VEL)); % lower bound
medmax = median(VEL)+(0.5*iqr(VEL)); % higher bound

% all to be binned 'MEDIUM' fall in between these values
% 'SLOW' Speed: lowest speed to medmin
% 'FAST' Speed: medmax to highest speed

speedBin = nan(size(VEL));

for i=1:length(VEL);
    if VEL(i) < medmin
        speedBin(i) = 1; % slow
    elseif VEL(i) > medmax
        speedBin(i) = 3; % fast
    else
        speedBin(i) = 2; % med
    end
end

VEL_isSlow = speedBin(:,1)==1;
VEL_isMed= speedBin(:,1)==2;
VEL_isFast= speedBin(:,1)==3;

slowData = data(VEL_isSlow,:);
medData = data(VEL_isMed,:);
fastData = data(VEL_isFast,:);


clear i medmax medmin speedBin

end

