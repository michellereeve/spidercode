function [smInterpData,der1,der2,spData] = SplineInterp_wSPAPS(data, time, sptol, spMode, newTime, plotMode, weights)
% [smInterpData, der1, der2, spData] = SplineInterp(data, time, sptol, spMode, newTime,plotMode,weights)
% Takes input data and interpolates any missing (NaN) values.
% Assumes missing values are indicated by 'NaN'
% Can be used on qualysis kinematic data, as long as the export option has
% been selected to use 'NaN' for missing values rather than zeros.
% Include optional input 'newTime' if the data is to be up-sampled as well
%
% Inputs
%   data: raw data to be smoothed and interpolated
%   time: time wave that correspondes to data (should have the same
%      number of points
%   sptol: smoothing factor which should be related to the
%       error variance of the input data, if known.
%       If this value is not known, use [] as a place holder if 'newTime' is an
%       input
%   spMode:  spline order: linear (1), cubic (2), quintic (3)
%       set to quintic spMode is not included as an input
%   newTime: Time wave to re-sample to, if desired. If this is not
%        included, the function will only fill in missing values in the original
%   time wave (denoted by 'NaN').
%   plotMode: if plotMode = 1, will plot results
%
% Outputs
%   interpData: smoothed and interpolated data
%   der1: first derivative of interpData
%   der2: second derivative of interpData
%   spData: spline function information


if ~exist('newTime', 'var') || isempty(newTime)
    newTime = time;
end

if ~exist('spMode', 'var') || isempty(spMode)
    spMode = 2;
end

if ~exist('sptol','var')  || isempty(sptol)
    sptol = 0;
end

if ~exist('plotMode','var')  || isempty(plotMode)
    plotMode = 0;
end

if ~exist('weights','var')  || isempty(weights)
    weightMode = 0;
else
    weightMode = 1;
end

nRows = length(newTime);
nCols = size(data,2);

smInterpData = NaN(nRows,nCols);
der1         = NaN(nRows,nCols);
der2         = NaN(nRows,nCols);

for i = 1:nCols
    t_data = data(:,i);
    t_time = time;   
    NaNpts = find(isnan(t_data));
    t_data(NaNpts) = [];
    t_time(NaNpts) = [];
    
    if weightMode == 0
        spData = spaps(t_time, t_data, sptol, spMode);
    else
        t_weights = weights;
        t_weights(NaNpts) = [];
        spData = spaps(t_time, t_data, sptol, spMode, t_weights);
    end
    
    smInterpData(:,i) = fnval(spData,newTime);
    der1(:,i) = fnval(fnder(spData,1),newTime);
    der2(:,i) = fnval(fnder(spData,2),newTime);
    
     smInterpData(NaNpts,i) = nan(length(NaNpts),1);
     der1(NaNpts,i) =  nan(length(NaNpts),1);
     der2(NaNpts,i) =  nan(length(NaNpts),1);

end

if plotMode == 1
    for i = 1:nCols
        figure
        plot(time,data(:,i),'r')
        hold on;
        plot(newTime,smInterpData(:,i),'b')
        title('Raw data (red), and Spline interpolated data (blue)')
    end

end

end





