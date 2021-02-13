%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function [measFilt, measFilt_minPhase] = getMeasurementFilter(sweep,invSweep,fs,filterSize,fsTarget)
%
% Function to get a filter which can be applied to synthesized full range
% BRIRs in order to adjust to the actual measurement range
%
% Output:
% measFilt              - Linear phase FIR filter in time domain
% measFilt_minPhase     - Minimum phase FIR filter in time domain
%
% Input:        
% sweep                 - Measurement sweep applied
% invSweep              - Inverse sweep applied for deconvolution
% filterSize            - Filter size in samples. Should be 2^n. Minimum 
%                         phase filter will have a length of filterSize/2
%                         Default: 2^12 taps
% fs                    - Sample rate of sweeps
% fsTarget              - Target sample rate of filter if
%                         up/downsampling should be applied. Leave empty if
%                         sample rate should not be changed.
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

function [measFilt, measFilt_minPhase] = getMeasurementFilter(sweep,invSweep,fs,filterSize,fsTarget)

if size(sweep,2) > size(sweep,1)
    sweep = sweep.';
end

if size(invSweep,2) > size(invSweep,1)
    invSweep = invSweep.';
end

if length(sweep) ~= length(invSweep)
    error('Sweep and inverse sweep must have the same length');
end

if nargin < 3
    error('Please specify sample rate of RIR');
end

if nargin < 4 || isempty(filterSize)
    filterSize = 2^12;
end

if nargin < 5 || isempty(fsTarget)
    fsTarget = fs; % no resampling;
end

%% Design filter

%Apply resampling if desired
if fs ~= fsTarget
    gComDiv = gcd(fs, fsTarget);
    p = fsTarget / gComDiv;
    q = fs / gComDiv;
    sweep = resample(sweep,double(p),double(q));
    invSweep = resample(invSweep,double(p),double(q));
end

%Get raw filter
L = length(sweep)+length(invSweep)-1;
NFFT = 2^nextpow2(L);
specSweep = fft(sweep,NFFT);
specInvSweep = fft(invSweep,NFFT);
measFilt = abs(specSweep .* specInvSweep);
measFilt = measFilt(1:end/2+1);

%Just assuming that the deconvolution is linear at 1kHz, we take the values
%at 1kHz octave range and set this to 0
fVec = linspace(0,fsTarget/2,length(measFilt));
[~,bin500] = min(abs(fVec-500));
[~,bin2000] = min(abs(fVec-2000));
measFilt = measFilt ./ (1-mean(measFilt(bin500:bin2000)));

%Resample to filterSize and smooth
measFilt = resample(measFilt,1,NFFT/filterSize);
measFilt = AKfractOctSmooth(measFilt,'amp',fsTarget,1);

%Build linear phase filter
measFilt = AKsingle2bothSidedSpectrum(measFilt);
measFilt = real(ifft(measFilt));
measFilt = circshift(measFilt, length(measFilt)/2-1);
measFilt = measFilt .* hann(length(measFilt));

%And get minimum phase version with rceps
[~,measFilt_minPhase] = rceps(measFilt);
measFilt_minPhase = measFilt_minPhase(1:end/2);
minPhaseWin = supdeq_win(length(measFilt_minPhase),[0,length(measFilt_minPhase)/2]);
measFilt_minPhase = measFilt_minPhase.*minPhaseWin;

end