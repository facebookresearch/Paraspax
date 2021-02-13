%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function psx = detectReflections(psx)
%
% Function to detect reflections in monaural IR according to [1].
%
% Output:
% psx                   - psx struct with list of detected reflections in
%                         psx.refDet
%
% Input:        
% psx                   - psx struct with required fields
%
% Dependencies: SUpDEq toolbox, AKtools, ITA toolbox
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

function psx = detectReflections(psx)

%Get various variables for further processing
fs = psx.fs;
rir = psx.rir;
tDirect = psx.par.mon.rirOnset_20dB; %TOA of direct sound in samples
directWin = round(psx.refDet.directWin*fs); %Direct sound window in samples
refWin = round(psx.refDet.refWin*fs); %Reflection window in samples
minTdist = psx.refDet.minTdist; %Minimum time distance between two reflections
prefPeaks = psx.refDet.prefPeaks; %Boolean if peaks should be prefered over energy for discarding
blockSize = round(psx.refDet.blockSize*fs); %Block size for reflection detection in samples
refFiltNFFT = psx.refDet.refFiltNFFT; %FFT size for reflection filters
refFiltLowFreq = psx.refDet.refFiltLowFreq; %Lowest detectable frequency for reflection filters (varies maximum window size)
mts  = round(psx.par.mon.mtAbel * 10^-3 * fs); %Mixing time according to abel in samples
mts2 = 2*mts; %2 times mixing time used for analyses;
percT = psx.refDet.percT; %Boolean if perceptional threshold should be applied for discarding reflections
percTrange = psx.refDet.percTrange; %Gain range of perceptual threshold in dB
nRef = psx.refDet.nRef; %Number of reflections to render / use for parametric description

%% Get RMS of direct sound

directTaps = tDirect-directWin(1):tDirect+directWin(2);
eDirect = rms(rir(directTaps));

%% Run through RIR up to mts2 and detect reflections (local energy 3 times higher than median in window)

numOfBlocks = ceil(mts2/blockSize)-1; %Number of blocks up to 2 * mixing time
eRIR = (rir).^2; %Squared rir

%Init arrays and blockIdx
peakVal = zeros(blockSize*numOfBlocks,1);
peakD = zeros(blockSize*numOfBlocks,1);
blockIdx = 0;
%Run through blocks
while(blockIdx <= numOfBlocks-1)
    %Get block of eRIR
    winSamples = blockIdx*blockSize+1:(blockIdx+1)*blockSize;
    block = eRIR(winSamples);
    medBlock = median(block);
    %Go through block and check if local energy is 3 times higher than
    %median, if yes mark sample
    for kk = winSamples(1):winSamples(end)
        if eRIR(kk) > 3*medBlock
            peakD(kk) = 1;
            peakVal(kk) = eRIR(kk);
        end
    end
    %Find max in peakVal window and mark as 100 (random id)
    [~,maxPeakId] = max(abs(peakVal(winSamples)));
    peakD(winSamples(1)-1+maxPeakId) = 100;
    %Iterate over blockIdx
    blockIdx = blockIdx+1;
end
%Get peak samples
peakS = find(peakD==100);

%Get rid of values in direct sound window
peakS(peakS <= directTaps(end) ) = [];

%% Sort and discard reflections

%Get energy of detected reflections in pre-defined reflection windows, peak
%level, and make first reflection list
for kk = 1:length(peakS)
    eRef(kk) = rms(rir(peakS(kk)-refWin(1):peakS(kk)+refWin(2)));
    pRef(kk) = abs(rir(peakS(kk)));
end
refListDet = [peakS,eRef',pRef'];

%Erase reflections which arive within minTdist after reflection with higher
%energy
dRef = diff(refListDet(:,1));
dRef = find(dRef <= minTdist*fs);
discardIDs = nan;
for kk = 1:length(dRef)

    %Energy based
    if ~prefPeaks
        [~,maxId] = max([refListDet(dRef(kk),2), refListDet(dRef(kk)+1,2)]);
        if maxId == 1
            discardIDs(kk) = dRef(kk)+1; %Discard the other one....
        else
            discardIDs(kk) = dRef(kk);
        end
    end
    
    %Peak based
    if prefPeaks
        [~,maxId] = max([refListDet(dRef(kk),3), refListDet(dRef(kk)+1,3)]);
        if maxId == 1
            discardIDs(kk) = dRef(kk)+1; %Discard the other one....
        else
            discardIDs(kk) = dRef(kk);
        end
    end
    
end
if ~isnan(discardIDs)
    refListDet(discardIDs,:) = [];
end


if ~prefPeaks %Sort according to energy level
    [~,idx] = sort(refListDet(:,2),'descend');
    refListDet = refListDet(idx,:);
else %Sort according to peak level
    [~,idx] = sort(refListDet(:,3),'descend');
    refListDet = refListDet(idx,:);
end

%% Compare amplitude of reflections against masking threshold - Can be applied optionally

%Get peak level or direct sound to adjust threshold curve respectively
peakLevelDirect = 20*log10(abs(max(abs(rir))));

%Absolute threshold, single LATERAL reflection with speech, IEC room
%Olive & Toole 1989, Fig 11
%Toole - Sound Reproduction - Fig 6.8/6.9 (Live IEC listening room)
l = [-13,  -17,  -12.5    -10,  -11,     -13,   -15,   -11,   -11.5,   -15  ]; %relative level
t = [0,    2.5,      5,   7.5,   10,      20,    30,    40,      60,    80  ]; %delay relative to direct sound in ms
if t(end) < t(end-1)
    t(end) = t(end-1) + 10;
end

%Converte and interpolate
l = l + peakLevelDirect;
t = round(t*fs/1000);
t = t+tDirect;
t = [1,t,mts2];
l = [l(1),l,l(end)];
if t(end) < t(end-1)
    t(end) = t(end-1) + 0.001*fs;
end
l = interp1(t,l,[1:t(end)]);
lLow = l-percTrange;

%Discard reflections below threshold lLow if desired
if percT
    
    %Determine peak levels of reflections. In Toole publications, the
    %author always write about "relative level". I assume here they talk
    %about peak level. The results would be a little different if compared
    %against the energy level, as there are sometimes various reflections
    %in the energy window not covered if checking only the peaks...
    discardIDs = nan(length(refListDet),1);
    for kk = 1:length(refListDet)
        peakRef = 20*log10(max(abs(rir(refListDet(kk,1)))));
        if peakRef < lLow(refListDet(kk,1))
            discardIDs(kk) = 1;
        end
    end
    
    %Erase reflections below threshold
    discardIDs = find(~isnan(discardIDs));
    refListDet(discardIDs,:) = [];
end

%Change nRef if bigger then list of reflections
if nRef > length(refListDet)
    nRef = length(refListDet);
end

%% Calculate weighting function to separate specular components from diffuse components in the early part of the RIR

wf = zeros(length(rir(1:mts2)),1);
for kk = 1:nRef
    win = refListDet(kk,1)-refWin(1):refListDet(kk,1)+refWin(2);
    wf(win) = 1;
end
%Set to 1 for direct sound
wf(1:tDirect+directWin(2)) = 1;

%Smooth out zeros accorindg to energy of RIR (no hard transitions to 0)
%Get smoothed RIR
sm = round(0.003*fs);
if mod(sm,2)
    sm = sm+1;
end
rirSM = convFFT(hann(sm),abs(rir(1:mts2)));
rirSM = rirSM(sm/2:length(rir(1:mts2))+sm/2-1);
rirSM = (rirSM/median(rirSM)-1); %-1 because of direct sound
rirSM = rescale(rirSM);

%Combine smoothed abs(rir) curve and weighting function to fill space
%between reflections
zeroValues = find(wf==0);
wf(zeroValues) = rirSM(zeroValues);
%Set to zero after last reflection
wf(max(refListDet(1:nRef,1))+refWin(2)+1:end) = 0;

%Smooth, shift, and normalize wf
sm = blockSize;
if mod(sm,2)
    sm = sm+1;
end
wf = convFFT(hann(sm),wf);
%And shift again
wf = wf(sm/2:length(rir(1:mts2))+sm/2-1);
%Normalize again 
wf = wf./max(abs(wf));
%Set to 1 for direct sound
wf(1:tDirect) = 1;

%% Get reflection filters

nOctChunks = 24;
nOctFilter = 1;

%%%% Direct sound
%Length for direct sound
lengthW1 = directWin(1)+directWin(2); %Length of selected direct sound window
fLw1 = round(fs/lengthW1);

lengthW3 = round(fs/refFiltLowFreq); %Length of longest window according to refFiltLowFreq
fLw3 = round(fs/lengthW3);

fLw2 = round(mean([fLw1,fLw3])); %Frequency in between both corner frequencies
lengthW2 = round(fs/fLw2); %Length of window in between

%Get magnitude spectrum of direct sound for different window sizes
dsw1 = rir(tDirect-directWin(1)+1:tDirect+directWin(2));
w1 = supdeq_win(length(dsw1),[floor(directWin(1)*1/6),floor(directWin(2)*1/6)]);
dsw1 = dsw1.*w1;
fdsw1 = abs(fft(dsw1,refFiltNFFT)/length(dsw1));
fdsw1 = fdsw1(1:end/2+1);
fdsw1 = AKfractOctSmooth(fdsw1,'amp',fs,nOctChunks);
fdsw1 = fdsw1./mean(fdsw1);

dsw2 = rir(tDirect-directWin(1)+1:tDirect+lengthW2-directWin(1));
w2 = supdeq_win(length(dsw2),[floor(directWin(1)*1/6),floor(lengthW2-directWin(1)*1/6)]);
dsw2 = dsw2.*w2;
fdsw2 = abs(fft(dsw2,refFiltNFFT)/length(dsw2));
fdsw2 = fdsw2(1:end/2+1);
fdsw2 = AKfractOctSmooth(fdsw2,'amp',fs,nOctChunks);
fdsw2 = fdsw2./mean(fdsw2);

dsw3 = rir(tDirect-directWin(1)+1:tDirect+lengthW3-directWin(1));
w3 = supdeq_win(length(dsw3),[floor(directWin(1)*1/6),floor(lengthW3-directWin(1)*1/6)]);
dsw3 = dsw3.*w3;
fdsw3 = abs(fft(dsw3,refFiltNFFT)/length(dsw3));
fdsw3 = fdsw3(1:end/2+1);
fdsw3 = AKfractOctSmooth(fdsw3,'amp',fs,nOctChunks);
fdsw3 = fdsw3./mean(fdsw3);

%Combine windows
fVecDirect = linspace(0,fs/2,refFiltNFFT/2+1);
[~,fLw1] = min(abs(fLw1-fVecDirect)); 
[~,fLw2] = min(abs(fLw2-fVecDirect)); 
[~,fLw3] = min(abs(fLw3-fVecDirect)); 

fdcomb = [fdsw3(1:fLw2);fdsw2(fLw2+1:fLw1);fdsw1(fLw1+1:end)];
fdcomb = AKfractOctSmooth(fdcomb,'amp',fs,nOctFilter);

%Normalize
fdcomb = fdcomb*(1/max(abs(fdcomb)));


%%%% Early reflections
%Length for early reflections
lengthW1 = refWin(1)+refWin(2); %Length of selected direct sound window
fLw1 = round(fs/lengthW1);

lengthW3 = round(fs/refFiltLowFreq); %Length of longest window according to refFiltLowFreq
fLw3 = round(fs/lengthW3);

fLw2 = round(mean([fLw1,fLw3])); %Frequency in between both corner frequencies
lengthW2 = round(fs/fLw2); %Length of window in between

%Get magnitude spectrum of reflections for different window sizes
for kk = 1:nRef
    refw1(:,kk) = rir(refListDet(kk)-refWin(1)+1:refListDet(kk)+refWin(2));
    refw2(:,kk) = rir(refListDet(kk)-refWin(1)+1:refListDet(kk)+lengthW2-refWin(1));
    refw3(:,kk) = rir(refListDet(kk)-refWin(1)+1:refListDet(kk)+lengthW3-refWin(1));
end
w1 = supdeq_win(length(refw1),[floor(refWin(1)*1/6),floor(refWin(2)*1/6)]);
refw1 = refw1.*w1;
frefw1 = abs(fft(refw1,refFiltNFFT)/length(refw1));
frefw1 = frefw1(1:end/2+1,:);
for kk = 1:nRef
    frefw1(:,kk) = AKfractOctSmooth(frefw1(:,kk),'amp',fs,nOctChunks);
end
frefw1 = frefw1./mean(frefw1);

w2 = supdeq_win(length(refw2),[floor(refWin(1)*1/6),floor(lengthW2-refWin(1)*1/6)]);
refw2 = refw2.*w2;
frefw2 = abs(fft(refw2,refFiltNFFT)/length(refw2));
frefw2 = frefw2(1:end/2+1,:);
for kk = 1:nRef
    frefw2(:,kk) = AKfractOctSmooth(frefw2(:,kk),'amp',fs,nOctChunks);
end
frefw2 = frefw2./mean(frefw2);

w3 = supdeq_win(length(refw3),[floor(refWin(1)*1/6),floor(lengthW3-refWin(1)*1/6)]);
refw3 = refw3.*w3;
frefw3 = abs(fft(refw3,refFiltNFFT)/length(refw3));
frefw3 = frefw3(1:end/2+1,:);
for kk = 1:nRef
    frefw3(:,kk) = AKfractOctSmooth(frefw3(:,kk),'amp',fs,nOctChunks);
end
frefw3 = frefw3./mean(frefw3);

%Combine windows
[~,fLw1] = min(abs(fLw1-fVecDirect)); 
[~,fLw2] = min(abs(fLw2-fVecDirect)); 
[~,fLw3] = min(abs(fLw3-fVecDirect)); 

frefcomb = [frefw3(1:fLw2,:);frefw2(fLw2+1:fLw1,:);frefw1(fLw1+1:end,:)];
%Normalize
frefcomb = frefcomb.*(1./max(frefcomb));

%Get difference spectra between direct and reflected sound in order to get
%reflection filters (Spectral difference between direct sound and reflected
%sound)
fDiff = frefcomb;
[~,highFreqBin] = min(abs(fVecDirect-10000));
for kk = 1:nRef
    fDiff(1:highFreqBin,kk) = frefcomb(1:highFreqBin,kk)./fdcomb(1:highFreqBin);
end

fDiff = fDiff./mean(fDiff);
for kk = 1:nRef
    fDiff(:,kk) = AKfractOctSmooth(fDiff(:,kk),'amp',fs,nOctFilter);
end
%Normalize
fDiff = fDiff.*(1./max(fDiff));

%% Write results in struct

psx.refDet.tDirect = tDirect; 
psx.refDet.eDirect = eDirect; %RMS of direct sound
psx.refDet.refListDet = refListDet; %Detected reflections
psx.refDet.wf = wf; %Weighting function
psx.refDet.dsFilt = fdcomb; %Direct sound filter
psx.refDet.refFilt = fDiff; %Reflection filter
psx.refDet.fVecRefFilt = fVecDirect; %Frequency vector
psx.refDet.percTcurve = l;
psx.refDet.percTcurveLow = lLow;

end