% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% script diffuseFieldEQ_KEMAR_Meas_ALFE_NoIC
%
% Script to create diffuse field response of measured HRTF set, diffuse 
% field filter, and binaural noise with diffuse field response of
% dummy head
%   
% Dependencies: SUpDEq toolbox, AKtools
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


%%
clear all;
close all;
clc;

load('KEMAR_HRIRs_Meas_ALFE_sfd_N35');

%Get sample rate
fs = HRIRs_sfd_N35.f(end)*2;

%Get HRTFs on a really dense grid (including weights)
sg = supdeq_lebedev(5810);

%Get HRIRs (not HRTFs as we want a higher resolution for the filter)
[hrirL,hrirR] = supdeq_getArbHRIR(HRIRs_sfd_N35,sg,'DEG',2,'ak');

%Define filter for diffuse field and diffuse field compensation filter
filterLength = 512;

%Get magnitude
hrtfL = abs(fft(hrirL,filterLength));
hrtfR = abs(fft(hrirR,filterLength));
hrtfL = hrtfL(1:end/2+1,:);
hrtfR = hrtfR(1:end/2+1,:);

%Weights are part of the grid
weights = sg(:,3);

%Get diffuse field response of dummy head
dfrL = sqrt(sum((hrtfL.^2 .* weights'),2));
dfrR = sqrt(sum((hrtfR.^2 .* weights'),2));
dfr  = mean([dfrL,dfrR],2);

%Set to local value above specified frequency kHz %Could be done nicer, but this usually
%works at such high frequencies without too much ringing in the filter
fMin = 200;
fMax = 16000;
fVec = linspace(0,fs/2,filterLength/2+1);
[~,tap] = min(abs(fMin-fVec));
dfr(1:tap) = dfr(tap+1);
[~,tap] = min(abs(fMax-fVec));
dfr(tap+1:end) = dfr(tap);

%Set to 0 dB at 0Hz
offset = -20*log10(abs(dfr(1)));
dfr = dfr*10^(offset/20);

%Get inverse
dfrInv = 1./dfr;

%Smooth both magnitude spectra (1/3 octave)
dfr_sm = AKfractOctSmooth(dfr,'amp',fs,3);
dfrInv_sm = AKfractOctSmooth(dfrInv,'amp',fs,3);

%Transform to linear phase filters
dfr         = AKsingle2bothSidedSpectrum(dfr);
dfr_sm      = AKsingle2bothSidedSpectrum(dfr_sm);
dfrInv      = AKsingle2bothSidedSpectrum(dfrInv);
dfrInv_sm   = AKsingle2bothSidedSpectrum(dfrInv_sm);

dfr = ifft(dfr);
dfr = circshift(dfr,filterLength/2-1);
dfr = dfr .* hann(filterLength);
dfr = dfr * 1/(max(abs(dfr)));

dfr_sm = ifft(dfr_sm);
dfr_sm = circshift(dfr_sm,filterLength/2-1);
dfr_sm = dfr_sm .* hann(filterLength);
dfr_sm = dfr_sm * 1/(max(abs(dfr_sm)));

dfrInv = ifft(dfrInv);
dfrInv = circshift(dfrInv,filterLength/2-1);
dfrInv = dfrInv .* hann(filterLength);
dfrInv = dfrInv * 1/(max(abs(dfrInv)));

dfrInv_sm = ifft(dfrInv_sm);
dfrInv_sm = circshift(dfrInv_sm,filterLength/2-1);
dfrInv_sm = dfrInv_sm .* hann(filterLength);
dfrInv_sm = dfrInv_sm * 1/(max(abs(dfrInv_sm)));

%Plot to check
supdeq_plotIR(dfr_sm,dfrInv_sm)

%And get minimum phase filter using the cepstral function instead of the
%less convenient way with uncle Hilbert
[~,dfr_min] = rceps(dfr); 
[~,dfr_sm_min] = rceps(dfr_sm); 
[~,dfrInv_min] = rceps(dfrInv); 
[~,dfrInv_sm_min] = rceps(dfrInv_sm);

dfr_min = dfr_min(1:end/2);
dfr_sm_min = dfr_sm_min(1:end/2);
dfrInv_min = dfrInv_min(1:end/2);
dfrInv_sm_min = dfrInv_sm_min(1:end/2);

%Plot to check
close all;
supdeq_plotIR(dfr_sm_min,dfrInv_sm_min,[],[],8)

%Save all in struct
DFR.dummyHead = 'KEMAR_Meas_ALFE';
DFR.fs = fs;
DFR.dfr = dfr;
DFR.dfrInv = dfrInv;
DFR.dfr_sm = dfr_sm;
DFR.dfrInv_sm = dfrInv_sm;
DFR.dfr_min = dfr_min;
DFR.dfrInv_min = dfrInv_min;
DFR.dfr_sm_min = dfr_sm_min;
DFR.dfrInv_sm_min = dfrInv_sm_min;
save('DFR_KEMAR_Meas_ALFE','DFR');

%% Generate binaural noise which can be used for BRIR synthesis

noiseLength = 10*fs; %10 seconds
binNoise = randn(noiseLength,2);

%Apply dfr to get characterstics of dummy head
binNoise = convFFT(dfr_sm_min,binNoise);

%Normalize
binNoise = binNoise * 0.99/max(max(abs(binNoise)));

%Save
save('BinauralWhiteNoise_KEMAR_Meas_ALFE_Artificial_NoIC','binNoise');
