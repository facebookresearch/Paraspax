%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function smax = slidingMax(x,winSize)
%
% Function to calculate max in a sliding window
%
% Output:
% srms                  - Sliding max vector with same size as input
%                         signal. Step size is one sample with max of 
%                         window according to winSize. Signal gets nan
%                         padded at beginning/end.
%
% Input:        
% x                     - Input signal
% winSize               - Window size in samples
%                         Default: 256 samples
%
% Dependencies: -
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.
%
% Copyright (c) Facebook, Inc. and its affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

function smax = slidingMax(x,winSize)

if nargin < 2 || isempty(winSize)
    winSize = 256;
end

if size(x,2) > size(x,1)
    x = x.';
end

%% Get max value in sliding window

% Zero pad first according to window size
if mod(winSize,2)
    winSize = winSize+1;
end

nanPad = nan(winSize/2,size(x,2));
x = [nanPad;x;nanPad];
sStart = winSize/2+1;
sEnd = length(x)-winSize/2;

smax = nan(size(x,1),size(x,2));
for kk = sStart:sEnd
    smax(kk,:) = max(abs((x(kk-winSize/2:kk+winSize/2-1,:))),[],'omitnan');
end
smax = smax(winSize/2+1:end-winSize/2,:);

end