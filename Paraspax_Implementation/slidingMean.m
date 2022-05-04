% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function smean = slidingMean(x,winSize)
%
% Function to calculate mean in a sliding window
%
% Output:
% srms                  - Sliding mean value vector with same size as input
%                         signal. Step size is one sample with mean of 
%                         window according to winSize. Signal gets nan
%                         padded at beginning/end leading.
%
% Input:        
% x                     - Input signal
% winSize               - Window size in samples
%                         Default: 32 samples
%
% Dependencies: -
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function smean = slidingMean(x,winSize)

if nargin < 2 || isempty(winSize)
    winSize = 32;
end

if size(x,2) > size(x,1)
    x = x.';
end

%% Get mean value in sliding window

% Zero pad first according to window size
if mod(winSize,2)
    winSize = winSize+1;
end

nanPad = nan(winSize/2,size(x,2));
x = [nanPad;x;nanPad];
sStart = winSize/2+1;
sEnd = length(x)-winSize/2;

smean = nan(size(x,1),size(x,2));
for kk = sStart:sEnd
    smean(kk,:) = mean(x(kk-winSize/2:kk+winSize/2-1,:),'omitnan');
end
smean = smean(winSize/2+1:end-winSize/2,:);

end