% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function srms = slidingRMS(x,winSize)
%
% Function to calculate rms in a sliding window
%
% Output:
% srms                  - Sliding rms value vector with same size as input
%                         signal. Step size is one sample with rms of 
%                         window according to winSize. Signal gets zero
%                         padded at beginning/end leading to slightly wrong
%                         rms values at beginning/end.
%
% Input:        
% x                     - Input signal
% winSize               - Window size in samples
%                         Default: 32 samples
%
% Dependencies: -
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function srms = slidingRMS(x,winSize)

if nargin < 2 || isempty(winSize)
    winSize = 32;
end

if size(x,2) > size(x,1)
    x = x.';
end

%% Get rms value in sliding window

% Zero pad first according to window size
if mod(winSize,2)
    winSize = winSize+1;
end

zeroPad = zeros(winSize/2,size(x,2));
x = [zeroPad;x;zeroPad];
sStart = winSize/2+1;
sEnd = length(x)-winSize/2;

srms = nan(size(x,1),size(x,2));
for kk = sStart:sEnd
    srms(kk,:) = rms(x(kk-winSize/2:kk+winSize/2-1,:));
end
srms = srms(winSize/2+1:end-winSize/2,:);

end