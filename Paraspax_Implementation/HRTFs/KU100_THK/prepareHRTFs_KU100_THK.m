%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% script prepareHRTFs_KU100_THK
%
% Example script to show how to prepare a SOFA file for usage in the
% Paraspax rendering. The script downloads a SOFA file, transform the HRTF set
% to SH domain, creates the diffuse field response of the measured HRTF set, 
% the diffuse field filter, and binaural noise with the diffuse field 
% response of the HRTF set.
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

%%
clear all;
close all;
clc;

%Download KU100 HRTF dataset from TH Köln
HRIRurl = 'https://zenodo.org/record/3928297/files/HRIR_L2702.sofa?download=1';
HRIRfolder = [pwd,'/'];
HRIRfilename = 'KU100_HRIR_L2702_THK.sofa';
disp('Downloading HRIR Dataset');
HRIRsave = websave([HRIRfolder HRIRfilename],HRIRurl);

%Load SOFA file and transform to SH domain
sofaObj = SOFAload(HRIRfilename);

%Transform to SH domain at sufficiently high order
%In this case we apply the least-squares solution with Tikhonov
%regularization (just in case). Pass a grid with weights to apply the
%closed-form solution. Also be aware about the SH order and whether your
%data are full-spherical. Data which are not full spherical should be
%transformed with Tikhonov regularization to reduce errors.

%N = 35 should be enough for this application. Transform with weights 
%(use supdeq_lebedev to get sampling grid) would lead to slightly different 
%results, especially for contralateral HRTFs. 
HRIRs_sfd = supdeq_sofa2sfd(sofaObj,35,[],2,'ak',10^-6); 

%Save HRIRs_sfd struct an use in 'demo.m' as HRTF file
save KU100_HRIRs_sfd_N35 HRIRs_sfd;

%% Get diffuse field response

%Get sample rate
fs = HRIRs_sfd.f(end)*2;

%Get HRTFs on a really dense grid (including weights)
sg = supdeq_lebedev(5810);

%Get HRIRs (not HRTFs as we want a higher resolution for the filter)
[hrirL,hrirR] = supdeq_getArbHRIR(HRIRs_sfd,sg,'DEG',2,'ak');

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
dfrL    = sqrt(sum((hrtfL.^2 .* weights'),2));
dfrR    = sqrt(sum((hrtfR.^2 .* weights'),2));
dfr   = mean([dfrL,dfrR],2);

%Set to 0 dB at 0Hz
offset = -20*log10(abs(dfr(1)));
dfr = dfr*10^(offset/20);
%Set to local value above specified frequency kHz %Could be done nicer, but this usually
%works at such high frequencies without too much ringing in the filter
fMax = 12000;
fVec = linspace(0,fs/2,filterLength/2+1);
[~,tap] = min(abs(fMax-fVec));
dfr(tap+1:end) = dfr(tap);

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
close all;
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
supdeq_plotIR(dfr_sm_min,dfrInv_sm_min)

%Save all in struct
DFR.dummyHead = 'KU100_THK';
DFR.fs = fs;
DFR.dfr = dfr;
DFR.dfrInv = dfrInv;
DFR.dfr_sm = dfr_sm;
DFR.dfrInv_sm = dfrInv_sm;
DFR.dfr_min = dfr_min;
DFR.dfrInv_min = dfrInv_min;
DFR.dfr_sm_min = dfr_sm_min;
DFR.dfrInv_sm_min = dfrInv_sm_min;
save('DFR_KU100_THK','DFR');

%% Generate binaural noise which can be used for BRIR synthesis

%Decide if binaural noise should be filterd with IC filter or not
applyICfilter = false;

noiseLength = 10*fs; %10 seconds
binNoise = randn(noiseLength,2);

if applyICfilter
    %Apply interaural coherence filter 
    binNoise = AKbinauralCoherence(binNoise,fs);
end

%Apply dfr to get characterstics of dummy head
binNoise = convFFT(dfr_sm_min,binNoise);

%Normalize
binNoise = binNoise * 0.99/max(max(abs(binNoise)));

%Save
if applyICfilter
    save('BinauralWhiteNoise_KU100_THK_Artificial','binNoise');
else
    save('BinauralWhiteNoise_KU100_THK_Artificial_NoIC','binNoise');
end
