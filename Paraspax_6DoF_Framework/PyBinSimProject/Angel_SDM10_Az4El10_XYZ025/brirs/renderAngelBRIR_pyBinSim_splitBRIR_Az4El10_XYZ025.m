%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% script renderAngelBRIR_pyBinSim_splitBRIR_Az4EL10_XYZ025
%
% Script to render BRIRs for plausibility study and demo application (for 
% more details, see [1]).
%
% Dependencies: Paraspax toolbox
%
% References:
% [1] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on 
% One Monaural Room Impulse Response,” Submitted for publication, 2020.
%
% Code written 2019/2020 by JMA, Johannes M. Arend.
%
% Copyright (c) Facebook, Inc. and its affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree. 

%Choose src and rec position
%Data in the git are only available for recPosID = 10 and srcPosID = 1:4
for recPosID = 10 %Position 10 used for demo and experiment
    for srcPosID = 1:4 %To render for all sources
        clear psx;

%% Positions 

%Loudspeaker-Positions 1-4
LSpositions = [6.5600    3.7550    1.0690;...
               9.1020    2.4190    1.4020;...
               10.6440   1.2710    1.8000;...
               6.7120    0.9940    1.7380];
%SDM-Positions 1-24
SDMpositions = [5.592  4.016  1.40;...
                5.592  3.016  1.40;...
                5.592  2.016  1.40;...
                5.592  1.016  1.40;...
                4.592  4.016  1.40;...
                4.592  3.016  1.40;...
                4.592  2.016  1.40;...
                4.592  1.016  1.40;...
                3.592  4.016  1.40;...
                3.592  3.016  1.40;... %SDM10 --> Used for demo and experiment
                3.592  2.016  1.40;... 
                3.592  1.016  1.40;...
                2.592  4.016  1.40;...
                2.592  3.016  1.40;...
                2.592  2.016  1.40;...
                2.592  1.016  1.40;...
                1.592  4.016  1.40;...
                1.592  3.016  1.40;...
                1.592  2.016  1.40;...
                1.592  1.016  1.40;...
                0.592  4.016  1.40;...
                0.592  3.016  1.40;...
                0.592  2.016  1.40;...
                0.592  1.016  1.40];
            
%Calculate srcDistance and srcDirection, required for correct direct sound
srcDistance = SDMpositions(recPosID,:) - LSpositions(srcPosID,:); srcDistance = sqrt(srcDistance(1)^2 + srcDistance(2)^2 + srcDistance(3)^2);

srcDirection = LSpositions(srcPosID,:) - SDMpositions(recPosID,:);
[srcDirection(1),srcDirection(2)] = cart2sph(srcDirection(1),srcDirection(2),srcDirection(3));
srcDirection(1) = 360-mod(360-rad2deg(srcDirection(1)),360);
srcDirection(2) = rad2deg(srcDirection(2));
srcDirection(3) = [];

%% Load RIR

if recPosID < 10
    fileName = ['2019-07-16_Angel_SDM_MeasLoc_0',num2str(recPosID),'_SourceNum_',num2str(srcPosID),'_az00_el00'];
else
    fileName = ['2019-07-16_Angel_SDM_MeasLoc_',num2str(recPosID),'_SourceNum_',num2str(srcPosID),'_az00_el00'];
end
load(fileName)

%Raw measurement
psx.raw.rirMeas = roomMeasure.roomIR;%roomMeasure.roomIRs(:,7);
psx.raw.fsMeas = roomMeasure.Fs;
psx.raw.sweep = roomMeasure.sweep;
psx.raw.invSweep = roomMeasure.invSweep;

clear roomMeasure;

%% Geometry

psx.geo.sourceDistance = srcDistance; %Source distance in m
psx.geo.sourceDirection = srcDirection; %Source DOA in degree
psx.geo.roomDimensions = true; %A flag to state if room dimensions are available or not

%If geometrical data ara available, insert here for image source model (ISM)
if psx.geo.roomDimensions
    psx.geo.roomLength = 11.733; %X Axis for ISM
    psx.geo.roomWidth = 4.741; %Y axis for ISM
    psx.geo.roomHeight = 4.619; %Z axis for ISM
    psx.geo.srcPos = LSpositions(srcPosID,:); %Source position. Origin at lower right corner for ISM
    psx.geo.recPos = SDMpositions(recPosID,:); %Receiver position. Origin at lower right corner for ISM
    psx.geo.ISMorder = 2; %Image source order
    psx.geo.prefN1 = false; %Boolean if 1st order image sources should be prefered over 2nd order
end

%Pseudo-randomized DOAs based on shoebox model which can be used if no ISM can be applied. 
%Reflection pattern will be rotated according to source position in spatialize function
if ~psx.geo.roomDimensions
    psx.geo.arbRefPattern = [0,-30; 319,-2; 10,-2; 85,-2; 342,4; 319,-53; 0,67; 274,-50; 266,-3; 58,-31; 145,0; 260,-3; 194,-3];
end

%% Postprocess rir, resample to new sample rate, and shift according to source distance

psx.fs  = 48000; %Sampling rate of entire processing
psx.c = 343; %Speed of sound in m/s

%Post-processing without normalization, so level differences between measurements 1:4 will remain
psx.rir = postProRIR(psx.raw.rirMeas,psx.raw.fsMeas,psx.fs,0,psx.geo.sourceDistance); 

%% Get various filters
%Filter can be applied for post-equalization of BRIRs or headphone-equalization 
%for auralizations.

%Headphone compensation
psx.filter.hpc = importdata('AKG_K1000_Closed_KEMAR.mat');
%Measurement filter that can be applied if the measurement signal was band limited
[psx.filter.measFilt,psx.filter.measFilt_min] = getMeasurementFilter(psx.raw.sweep,psx.raw.invSweep,psx.raw.fsMeas,4096,psx.fs);

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

%Used default of -20 dB --> No reverberation level estimation

%% Perform spatialization with ISM or with random directions

psx = spatialize(psx);

%% Synthesize BRIR based on monaural RIR and estimated parameters

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
psx.brir.sg = getEquiangleGrid(4,10,true,0,50);

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

%% Save psx

if recPosID < 10
    psxName = ['psx_Angel_RecPos_0',num2str(recPosID),'_SrcPos_',num2str(srcPosID)];
else
    psxName = ['psx_Angel_RecPos_',num2str(recPosID),'_SrcPos_',num2str(srcPosID)];
end

folderName = psxName;
allFolderNames{1} = folderName;
mkdir(folderName)

%Save BRIRs
attenuation = 10^(-50/20); %Global attenuation estimated in various trials
safetyMargin = round(1.5*psx.fs/1000); %1.5 ms safety margin
winIn = 32;
winOut = 32;
renderBlockSize = 512; %Blocksize in PyBinSim
startSample = psx.par.mon.rirOnset_20dB-safetyMargin;

%Post-process synthesized BRIR and split for PyBinSim
brir = psx.brir.synthBRIR_spat;
%Get rid of propagation delay
brir = brir(startSample:end,:,:);
brir = brir .* supdeq_win(length(brir),[winIn 0]);
%Cut and fade to multiple of blockSize
%Late BRIR will have different length dependent on used IR. But not 
%important here as just one late BRIR is gonne be used anyway. However, if
%multipe late BRIRs are used, this has to be globally adjusted (or the late
%BRIRs are zero-padded / truncated to the same size in a later step)
endSample = mod(length(brir),renderBlockSize); 
endSample = length(brir)-endSample;
brir = brir(1:endSample,:,:);
brir = brir .* supdeq_win(length(brir),[0 winOut]);

%Cut after last dynamic reflexion (being dividable by blockSize)
%mtCut = max(psx.spat.refListSpat.toa) + (psx.refDet.refWin(end)*psx.fs + (size(psx.brir.HRTFs.f,2)*2-1)/psx.brir.HRTFs.FFToversize);
%Cut after mixing time
%EDIT - NICE TO HAVE, BUT MT DIFFERENT FOR EACH SOURCE - NOT PRACTICAL HEAR
%EVALUATED ALL maxTOA and applied mtCUT based on that!
%mtCut = ceil(psx.par.mon.mtAbel * psx.fs / 1000);
%modVal = mod(mtCut,renderBlockSize);
%mtCut = mtCut-modVal+renderBlockSize;

%In the end, we defined mtCut based evaluation of all sources (512*7) 
%The value will be enough, especially as BRIRs get additionally shifted 
% to get rid of propagation delay
mtCut = 3584; 

%Split BRIRs
brirEarly = brir(1:mtCut,:,:);
brirLate = brir(mtCut+1:end,:,1);

%Save as wav file for each head orientation
for kk = 1:size(brirEarly,3)
    sg = psx.brir.sg;
    fileName = [folderName,'/az',num2str(sg(kk,1)),'el',num2str(sg(kk,2)),'.wav'];
    audiowrite(fileName,(squeeze(brirEarly(:,:,kk))*attenuation),psx.fs,'BitsPerSample',32);
end
fileNameRev = [folderName,'/brirLate.wav'];
audiowrite(fileNameRev,(brirLate*attenuation),psx.fs,'BitsPerSample',32);

%% Perform extrapolation, render BRIRs, and save in psx structs

%Position shifts of 0.25 m on the x and y axis starting from the origin (P10)
xShift = [-2.25:0.25:2.25]';
yShift = [-2.25:0.25:1.25]'; %0.25m more then grid
xgrid = repmat(xShift,length(yShift),1);
ygrid = repmat(yShift,1,length(xShift)); ygrid = ygrid'; ygrid = ygrid(:);
cartGrid = [xgrid,ygrid,zeros(length(xgrid),1)];

%% Save grids

save('cartGrid','cartGrid');

sphGrid = psx.brir.sg;
save('sphGrid','sphGrid');

%% Perform extrapolation - translational shift - to new position

for nShift = 1:length(cartGrid)
    
    psx.extrap.lps = cartGrid(nShift,:);
    
    %Apply parameter extrapolation
    psx = extrapolate(psx);
    
    %Synthesize BRIRs
    psx.brir.synthMode = 'extrap';
    psx = synthesizeBRIR(psx);
    
    %Can be used to save psx in new folder 
    %If needed, uncomment, but data amount is huge!
    psxNameExtrap = [psxName,'_X_',num2str(psx.extrap.lps(1)),'_',num2str(psx.extrap.lps(2)),'_',num2str(psx.extrap.lps(3))];
    folderName = psxNameExtrap;
    allFolderNames{nShift+1} = folderName;
    mkdir(folderName)
    %save([folderName,'/',psxNameExtrap],'psx');

    %Save BRIRs
    startSample = psx.extrap.refListSpat.tDirect-safetyMargin; %Get extrap tDirect here!
    brir = psx.brir.synthBRIR_extrap; %Get extrap BRIR here!
    %Get rid of propagation delay
    brir = brir(startSample:end,:,:);
    brir = brir .* supdeq_win(length(brir),[winIn 0]);
    %Cut and fade to multiple of blockSize
    brir = brir(1:endSample,:,:);
    brir = brir .* supdeq_win(length(brir),[0 winOut]);
    
    %Error check
    if mtCut <= max(psx.extrap.refListSpat.toa)
        error('Extrapolated reflection after mixing time');
    end
    
    %Split BRIRs
    brirEarly = brir(1:mtCut,:,:);
    brirLate = brir(mtCut+1:end,:,1);

    %Save as wav file for each head orientation
    for kk = 1:size(brirEarly,3)
        sg = psx.brir.sg;
        fileName = [folderName,'/az',num2str(sg(kk,1)),'el',num2str(sg(kk,2)),'.wav'];
        audiowrite(fileName,(squeeze(brirEarly(:,:,kk))*attenuation),psx.fs,'BitsPerSample',32);
    end
    
    fileNameRev = [folderName,'/brirLate.wav'];
    audiowrite(fileNameRev,(brirLate*attenuation),psx.fs,'BitsPerSample',32);
end

%% Write pybinsim filterList txt

allFolderNames = allFolderNames(:,2:end); %Dont use first non extrapolated brir and just add to 0 0 0 later
for kk = 1:length(allFolderNames) 
   wavPath{kk} = ['brirs/',allFolderNames{kk}];
end

az = sphGrid(:,1);
el = sphGrid(:,2); %Work with colatiude because of OSC 
xID = [0:length(wavPath)-1];
sourceID = srcPosID-1; %Starts at 0 in pybinsim

fid=fopen(['filterList_FRL_Angle_SDM10_SrcPos_',num2str(srcPosID),'.txt'],'w');
for nID = 1:length(xID)
    for nSG = 1:length(az)
            txtLine = ['FILTER ',num2str(az(nSG)),' ',num2str(el(nSG)),' ',num2str(xID(nID)),' ',num2str(sourceID),' 0 0 ',wavPath{nID},'/az',num2str(az(nSG)),'el',num2str(el(nSG)),'.wav\n'];
            fprintf(fid, txtLine);
    end
end

%Add late reverb
txtLine = ['LATEREVERB 0 0 0 0 0 0 latereverb/brirLate.wav\n'];
fprintf(fid, txtLine);
%Add hp filter
txtLine = ['HPFILTER hpirs/AKG_K1000_Closed_KEMAR_2Channel.wav'];
fprintf(fid, txtLine);

fclose(fid);

    end
end