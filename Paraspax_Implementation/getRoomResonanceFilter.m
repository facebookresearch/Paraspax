%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function [resFilt, resFilt_minPhase] = getRoomResonanceFilter(rir,filterSize,maxGain,meanRange,lowFreqRange,fs)
%
% Function to get a filter which can be applied to synthesized full range
% BRIRs in order to compensate for room resonances. Can be handy to flatten
% the lower frequency range if the monaural RIR was measured at a position
% with strong room modes.
%
% Output:
% resFilt               - Linear phase FIR filter in time domain
% resFilt_minPhase      - Minimum phase FIR filter in time domain
%
% Input:        
% rir                   - Single channel RIR
% filterSize            - Filter size in samples. Should be 2^n. Minimum 
%                         phase filter will have a length of filterSize/2
%                         Default: 2^12 taps
% maxGain               - Maximum increase/decrease in dB
%                         Default: 6 dB
% meanRange             - 2 element vector defining lower and upper corner 
%                         frequency or range used to estimate the mean level 
%                         of the magnitude response (in Hz)
%                         Default: [100 2000] (Hz)
% filterRange           - 2 element vector defining lower and upper corner
%                         frequency of range where EQing is applied (in Hz)
%                         Default: [100 500] (Hz)
% fs                    - Sample rate of RIR
%                         Default: 48000
%
% Dependencies: SUpDEq toolbox, AKtools
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.
%
% Copyright (c) Facebook, Inc. and its affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

function [resFilt, resFilt_minPhase] = getRoomResonanceFilter(rir,filterSize,maxGain,meanRange,filterRange,fs)

if size(rir,2) > size(rir,1)
    rir = rir.';
end

if nargin < 2 || isempty(filterSize)
    filterSize = 2^12;
end

if nargin < 3 || isempty(maxGain)
    maxGain = 6;
end

if nargin < 4 || isempty(meanRange)
    meanRange = [100 2000];
end

if nargin < 5 || isempty(filterRange)
    filterRange = [100 500];
end

if nargin < 6 || isempty(fs)
    fs = 48000;
end

%% Design filter

%Get smoothed magnitude of RIR first
NFFT  = 2^nextpow2(length(rir));
magSM = abs(fft(rir,NFFT));
magSM = magSM(1:end/2+1);
magSM = AKfractOctSmooth(magSM,'amp',fs,6);

%Define frequency vector
fVec = linspace(0,fs/2,length(magSM));

%Get bins for mean value
[~,lowBinMeanRange] = min(abs(fVec-meanRange(1)));
[~,highBinMeanRange] = min(abs(fVec-meanRange(2)));

%Get mean value
meanVal = mean(magSM(lowBinMeanRange:highBinMeanRange));

%Get bins for filter
[~,lowBinFilterRange] = min(abs(fVec-filterRange(1)));
[~,highBinFilterRange] = min(abs(fVec-filterRange(2)));

%Design filter 
magFilter = ones(length(fVec),1)*meanVal;
magFilter(lowBinFilterRange:highBinFilterRange) = magSM(lowBinFilterRange:highBinFilterRange);
shift = -20*log10(meanVal);
magFilter = magFilter*(10^(shift/20));
magFilter = 1./magFilter;

%Adjust according to max gain
magFilter(magFilter > 10^(maxGain/20)) = 10^(maxGain/20);
magFilter(magFilter < 10^(-maxGain/20)) = 10^(-maxGain/20);

%Smooth again as magnitude response might have some edges
magFilter = AKfractOctSmooth(magFilter,'amp',fs,6);

%Resample to filter size
magFilter = resample(magFilter,1,NFFT/filterSize);

%Design linear phase filter
resFilt = AKsingle2bothSidedSpectrum(magFilter);
resFilt = real(ifft(resFilt));
resFilt = circshift(resFilt,length(resFilt)/2-1);
resFilt = resFilt.*hann(length(resFilt));

%Design minimum phase version with rceps
[~,resFilt_minPhase] = rceps(resFilt);
resFilt_minPhase = resFilt_minPhase(1:end/2);
minPhaseWin = supdeq_win(length(resFilt_minPhase),[0,length(resFilt_minPhase)/2]);
resFilt_minPhase = resFilt_minPhase.*minPhaseWin;

end