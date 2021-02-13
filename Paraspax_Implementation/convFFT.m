%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function c = convFFT(s,f)
%
% Just a simple function to perform (2-channel) convolution in frequency domain.
%
% Output:
% c                     - Convolution result
%
% Input:        
% s                     - Audio signal to be filtered
% f                     - Single- or double-channel filter which should be
%                         applied to signal s
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

function c = convFFT(s,f)

if size(s,2) > size(s,1)
    s = s.';
    warning('Flipped dimensions of signal');
end

if size(f,2) > size(f,1)
    f = f.';
    warning('Flipped dimensions of filter');
end

L = length(s)+length(f)-1;
NFFT = 2^nextpow2(L);

%For two channel filter
if size(f,2) > 1
    s = [s,s];
end

c = ifft(fft(s,NFFT) .* fft(f,NFFT));
c = c(1:L,:);