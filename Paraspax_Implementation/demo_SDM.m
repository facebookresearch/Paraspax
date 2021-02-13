%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% script demo_SDM
%
% Script to run a demo of the Paraspax toolbox [1] for synthesizing a BRIR 
% based on a monaural RIR and a DOA vector estimated with the SDM method [2]
%
% Dependencies: Paraspax toolbox, SUpDEq toolbox, AKtools, ITA Toolbox, Resources
%
% References:
% [1] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on 
% One Monaural Room Impulse Response,” Submitted for publication, 2020.
%
%﻿[2] S. V. Amengual Garí, J. M. Arend, P. Calamia, and P. W. Robinson, 
% “Optimizations of the Spatial Decomposition Method for 
% Binaural Reproduction (in press),” J. Audio Eng. Soc., 2020.
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.
%
% Copyright (c) Facebook, Inc. and its affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

%% Download required toolboxes and place in thirdParty folder

%SUpDEq toolbox: https://github.com/AudioGroupCologne/SUpDEq
%AKtools: https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/
%ITA toolbox: https://git.rwth-aachen.de/ita/toolbox

%% Load SRIR measurements

clear all;
close all;
clc;

%Add paths
addpath(genpath(pwd))

%Load SDM measurement struct
load('SRIR_data')

%% Postprocess RIR, normalize and shift according to source distance

%Source distance must always be available, independent of the
%spatialization mode. The source distance could also be encoded from the
%measurement (if appropriately latency-compensated), but in this case we
%just use data we captured with the OptiTrack system.

psx.geo.sourceDistance = 5.5422; %Source distance in m

psx.fs = 48000; %Sampling rate of entire processing
psx.c = 343; %Speed of sound in m/s
psx.geo.doa = SRIR_data.DOA; %Smoothed DOA vector
psx.geo.roomDimensions = false; %No room dimensions available. Spatialize function checks this flag

%Normal post-processing cannot be applied here as DOA vector is related to RIR%
%Only normalization and shift to sourceDistance applied, including shift of
%the DOA vector
[psx.rir,psx.geo.doa] = postProRIR(SRIR_data.P_RIR,psx.fs,[],1,psx.geo.sourceDistance,psx.geo.doa);

%% Get various filters
%Filter can be applied for post-equalization of BRIRs or headphone-equalization 
%for auralizations.

%Headphone compensation - Use other HPCF if other HRTF set is used!
psx.filter.hpc = importdata('Sennheiser-HD600_LessBass_KEMAR.mat');
%Room resonance filter that can be applied to compensate for resonances
[psx.filter.resFilt,psx.filter.resFilt_min] = getRoomResonanceFilter(psx.rir,4096,6); %Maximum boost set to 6dB

%% Get basic monaural room acoustic parameters in octave bands and broad band

psx.par.mon = getMonPar(psx.rir,psx.fs,1);

%% Define parameters for reflection detection and run detection function

%Window (in s) around detected direct sound
psx.refDet.directWin = [0.0005, 0.001];
%Window (in s) around detected reflection
psx.refDet.refWin = [0.0005, 0.001];
%Minimal distance (in s) between reflection peaks. Higher peak will be chosen if distance is smaller
psx.refDet.minTdist = 0.00075;
%Boolean if prefer peaks over energy in minTdist window
psx.refDet.prefPeaks = true;
%Block size (in s) for reflection detection
psx.refDet.blockSize = 0.001;
%FFT size for reflection filters
psx.refDet.refFiltNFFT = 512;
%Lowest detectable frequency for reflection filters (varies maximum window size)
psx.refDet.refFiltLowFreq = 50;
%Number of loudest reflections to be rendered / parameterized
psx.refDet.nRef = 10;
%Boolean to apply comparison against perceptional threshold according to
%Olive&Toole(1998)
psx.refDet.percT = false;
%Gain range of perceptual threshold in dB
psx.refDet.percTrange = 10;

%Run function to detect reflections and write results in psx struct
psx = detectReflections(psx);

%% Estimate level of reverberation

%Define method to estimate reverberation level.
%EDC - EDC curve will be transformed to level curve and shifted according
%to level of direct sound. Values will be estimated based on EDC curve.
%RMS - Sliding RMS window will be applied to estimate RMS level curve.
%Values will be estimated based on RMS curve.
%MAX - Sliding MAX window will be applied to estimate MAX level curve.
%Values will be estimated based on MAX curve (as described in [1]).
%'EDC_RMS_MAX' - All methods are gonne be applied and can be compared
% In this case, the BRIR synthesis favors the 'MAX' method, as this
% provided the best results (see [1])
psx.rev.method = 'EDC_RMS_MAX'; 
%The results are saved in psx.rev with various possible reverberation
%curves. Check by plotting later.
psx = estimateReverbLevel(psx);

%% Perform spatialization using DOAs estimated with SDM

psx = spatialize(psx);

%% Perform extrapolation - translational shift of listener to new position

%X Y Z for listener translation in meter
%Example: Move 2 m on X axis to the back, 1 m on Y axis to the
%right (origin of coordinate system in lower right corner).
%Results in a translational shift from Pos. 10 to Pos. 19 (cf., [1])
psx.extrap.lps = [-2,-1,0]; 

%Apply extrapolation and rotation according to head orientation
psx = extrapolate(psx);

%% Synthesize BRIR based on monaural RIR and estimated parameters

%%%% HRTFs %%%%
%Change the following 3 lines of code if you want to use a different HRTF set. 
%The Paraspax synthesis requires the HRTF set in SH domain, eventually the
%respective DFR, and the binaural white noise shaped with the diffuse field
%response of the respective HRTF set. The script 'prepareHRTFs_KU100_THK'
%in the folder 'HRTFs/KU100_THK' provides step by step instructions on how
%to convert an HRTF set in SOFA format to SH domain and on how to generate
%the diffuse field response and the binaural white noise
%%%%%%%%%%%%%%%

%HRTFs in SH domain
psx.brir.HRTFs = importdata('KEMAR_HRIRs_MeasSS2_ALFE_sfd_N35.mat');

%Load corresponding diffuse field response / common transfer function
psx.brir.DFR = importdata('DFR_KEMAR_MeasSS2_ALFE.mat');

%Load binaural white noise for synthesis of binaural diffuse reverberation
%Should match the HRTF dataset. Here without IC filter as it is gonne be
%applied in a later processing step.
psx.brir.binNoise = importdata('BinauralWhiteNoise_KEMAR_MeasSS2_ALFE_Artificial_NoIC.mat');

%Define blockSize of rir for convolution in samples
psx.brir.binRevBlockSizeRir = 32; %Provides the best results, based on informal comparisons to reference BRIRs
%Define blockSize of white noise for convolution in samples
psx.brir.binRevBlockSizeNoise = 128; %128 or 256 provides best results, based on informal comparisons to reference BRIRs

%Load source directivity if available
psx.brir.sourceDirectivity = importdata('Genelec_8020_directivity.mat');
%Define filter length of directivity filters
psx.brir.directivityFilterLength = 256;

%Define if measurement filter should be applied or not
psx.brir.applyMeasFilt = false; %Default - false

%Define if diffuse reverb should be filtered with a highpass (7th order
%butterworth at 200 Hz)
psx.brir.hpBinRev = false; %Default - false

%Define if Butterworth xover should be used and low frequency component of
%monaural RIR for the final BRIR
psx.brir.applyXover = false; %Default - false

%Define if interaural coherence filter (IC) should be applied afterwards to
%all BRIRs. Makes sense in combination with non-IC-filtered diffuse noise
psx.brir.applyIC = true; %Default - true

%Set order N and corner frequency of butterworth highpass/lowcut. If set to
%0, no filter will be applied. Default order N = 3
%Applied here as RIR measurements are band limited to 0.1-20 kHz
psx.brir.butterFc = 100;
psx.brir.butterN = 4;

%Define if DRR adjustment (energy match after spatialization of direct
%sound and early reflections) should be rms based or based on magnitude
%between 200Hz and 15kHz. String 'rms' or 'mag'
psx.brir.drrMatch = 'mag'; %Default mag

%Define sampling grid for BRIRs
%Have a look at the sampling grid with supdeq_plotGrid
psx.brir.sg = getEquiangleGrid(30,90,true,0,30);

%Define coordinate system of sampling grid: 
%global (=relative to sound sources) or local (=head-related)
%Global will result in sg(0,0) facing the sound source (relative angles)
%Local will result in sg(0,0) looking straight forward and the direct sound
%impinging from whatever direction it comes in the first place
psx.brir.sgSystem = 'global';

%Define if measured / spatialized or extrapolated parameters should be
%synthesized. String with 'spat' or 'extrap'
psx.brir.synthMode = 'spat';
psx = synthesizeBRIR(psx); %Result in psx.brir.synthBRIR_spat

%Also synthesize for the example extrapolated position
psx.brir.synthMode = 'extrap';
psx = synthesizeBRIR(psx); %Result in psx.brir.synthBRIR_extrap

%% Play example

git = audioread('Flamenco2_U89.wav');
%Apply HPeq
git = convFFT(git,psx.filter.hpc.minPhase);

brir_spat_sp2 = squeeze(psx.brir.synthBRIR_spat(:,:,2));
brir_extrap_sp2 = squeeze(psx.brir.synthBRIR_extrap(:,:,2));

git_spat_sp2 = convFFT(git,brir_spat_sp2);
git_extrap_sp2 = convFFT(git,brir_extrap_sp2);

sound(git_spat_sp2,psx.fs);
%sound(git_extrap_sp2,psx.fs);