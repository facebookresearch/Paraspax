% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function [combFilter, combFilter_minPhase] = combineFilters(filter1,filter2,filterSize)
%
% Function to combine 2 filters to one filter kernel. Only magnitude
% spectrum of filters gets combined, not the phase!
%
% Output:
% combFilter            - Combined linear phase FIR filter in time domain
% combFilter_minPhase   - Combined minimum phase FIR filter in time domain
%
% Input:        
% filter1               - Filter 1 in time domain, ideally linear phase
% filter2               - Filter 2 in time domain, ideally linear phase
% filterSize            - Filter size in samples. Should be 2^n. Minimum 
%                         phase filter will have a length of filterSize/2
%                         Default: 2^12 taps
%
% Dependencies: SUpDEq toolbox, AKtools
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function [combFilter, combFilter_minPhase] = combineFilters(filter1,filter2,filterSize)

if size(filter1,2) > size(filter1,1)
    filter1 = filter1.';
end

if size(filter2,2) > size(filter2,1)
    filter2 = filter2.';
end

if nargin < 3 || isempty(filterSize)
    filterSize = 2^12;
end

%% Design combined filter

%Get magnitude response of combined filter
magFilter = convFFT(filter1,filter2);
NFFT = 2^nextpow2(length(magFilter));
magFilter = abs(fft(magFilter,NFFT));
magFilter = magFilter(1:end/2+1);

%Resample to filter size
magFilter = resample(magFilter,1,NFFT/filterSize);

%Design linear phase filter
combFilter = AKsingle2bothSidedSpectrum(magFilter);
combFilter = real(ifft(combFilter));
combFilter = circshift(combFilter,length(combFilter)/2-1);
combFilter = combFilter.*hann(length(combFilter));

%Design minimum phase version with rceps
[~,combFilter_minPhase] = rceps(combFilter);
combFilter_minPhase = combFilter_minPhase(1:end/2);
minPhaseWin = supdeq_win(length(combFilter_minPhase),[0,length(combFilter_minPhase)/2]);
combFilter_minPhase = combFilter_minPhase.*minPhaseWin;

end