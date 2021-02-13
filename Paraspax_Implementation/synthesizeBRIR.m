%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function psx = synthesizeBRIR(psx)
%
% Function to synthesize BRIRs based on the parameterized monaural RIR.
% BRIRs can be synthesized for measurement point or extrapolated point.
% Passing a array with head orientations synthesized BRIRs for all passed
% head orientations.
%
% Output:
% psx                   - psx struct with synthesized BRIRs in psx.brir field
%
% Input:        
% psx                   - psx struct with required fields in psx.brir
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

function psx = synthesizeBRIR(psx)

%Get various variables for further processing
fs  = psx.fs;
rir = psx.rir;
nRef = psx.refDet.nRef;
mts  = round(psx.par.mon.mtAbel * 10^-3 * fs);
mts2 = 2*mts; 
wf = psx.refDet.wf;
hrtfs = psx.brir.HRTFs;
sgSH = psx.brir.sg; %SH coordinates not really used anymore
sg = sgSH; sg(:,2) = 90-sg(:,2); %Transform directly
nHeadOr = size(sgSH,1);
synthMode = psx.brir.synthMode;

if isfield(psx.brir,'sgSystem')
    sgSystem = psx.brir.sgSystem;
else
    psx.brir.sgSystem = 'global';
    sgSystem = 'global';
end
if isfield(psx.brir,'sourceDirectivity')
    sourceDirectivity = psx.brir.sourceDirectivity;
    directivityFilterLength = psx.brir.directivityFilterLength;
    applyDirectivity = true;
end
if isfield(psx.filter,'measFilt')
    measFilt = psx.filter.measFilt;
    measFiltMin = psx.filter.measFilt_min;
end
if isfield(psx.brir,'drrMatch')
    drrMatch = psx.brir.drrMatch;
else
    psx.brir.drrMatch = 'mag';
    drrMatch = 'mag';
end
if isfield(psx.brir,'applyMeasFilt')
    applyMeasFilt = psx.brir.applyMeasFilt;
else
    psx.brir.applyMeasFilt = false;
    applyMeasFilt = false;
end
if isfield(psx.brir,'hpBinRev')
    hpBinRev = psx.brir.hpBinRev;
else
    psx.brir.hpBinRev = false;
    hpBinRev = false;
end
if isfield(psx.brir,'applyXover')
    applyXover = psx.brir.applyXover;
else
    psx.brir.applyXover = false;
    applyXover = false;
end
if isfield(psx.brir,'applyIC')
    applyIC = psx.brir.applyIC;
else
    psx.brir.applyIC = true;
    applyIC = true;
end
if isfield(psx.brir,'butterFc')
    butterFc = psx.brir.butterFc;
    if isfield(psx.brir,'butterN')
        butterN = psx.brir.butterN;
    else
        butterN = 3;
    end
else
    psx.brir.butterFc = 0;
    psx.brir.butterN = 0;
    butterFc = 0;
    butterN = 0;
end
if isfield(psx,'rev')
    %Get value according to applied estimation method
    %If all methods were applied, use 'MAX' method as it seems to provide the best estimates
    if strcmp(psx.rev.method,'EDC')
        revAmplitude = psx.rev.edcMethod.interpLine_revAmp_firstRef;
    elseif strcmp(psx.rev.method,'RMS')
        revAmplitude = psx.rev.rmsMethod.interpLine_revAmp_firstRef;
    elseif strcmp(psx.rev.method,'MAX')
        revAmplitude = psx.rev.maxMethod.interpLine_revAmp_firstRef;
    elseif strcmp(psx.rev.method,'EDC_RMS_MAX')
        revAmplitude = psx.rev.maxMethod.interpLine_revAmp_firstRef;
    end
else
    revAmplitude = 0.1; %-20dB as default
end
psx.brir.revAmplitude = revAmplitude;

dfr = psx.brir.DFR;
toaDirect = psx.spat.refListSpat.tDirect;
directSE = psx.spat.refListSpat.SE_Direct;
azelDirect = psx.spat.refListSpat.AzEl_Direct;
refSE = psx.spat.refListSpat.selectSE;
azelRef = psx.spat.refListSpat.selectAzEl;
srcAzElDirect = psx.spat.refListSpat.SrcAzEl_Direct;
eDirect = psx.refDet.eDirect;
eRef = psx.refDet.refListDet(1:nRef,3);

if strcmp(synthMode,'extrap')
    toaDirect_ex = psx.extrap.refListSpat.tDirect;
    directSE_ex = psx.extrap.refListSpat.SE_Direct;
    azelDirect_ex = psx.extrap.refListSpat.AzEl_Direct;
    ampFactorDirect_ex = psx.extrap.refListSpat.ampFactorDirect;
    refSE_ex = psx.extrap.refListSpat.selectSE;
    azelRef_ex = psx.extrap.refListSpat.selectAzEl;
    refAmpFactor_ex = psx.extrap.refListSpat.ampFactor;
    srcAzElDirect_ex = psx.extrap.refListSpat.SrcAzEl_Direct;
    if ~isfield(psx.brir,'brirSpecGain')
        error('Global gain value for specular part of BRIR required. Please run spat synthesis mode first');
    else
        brirSpecGain = psx.brir.brirSpecGain;
    end
end

%% Filter RIR

if applyXover
    %Interaural Correlation below 200Hz almost one --> apply monaural RIR below
    %200 Hz.
    fCut = 200;

    %Get Linkwitz Riley of 7th order
    [bLowpass,aLowpass] = butter(7,fCut/fs*2,'low'); 
    [bHighpass,aHighpass] = butter(7,fCut/fs*2,'high'); 
    
    %Get different parts of RIR
    rirHigh = filter(bHighpass,aHighpass,rir);
    rirLow  = filter(bLowpass,aLowpass,rir);
else 
   rirHigh = rir; 
end

%% Generate binaural reverberation based on rir

%Check if already available and generate only if required
if ~isfield(psx.brir,'binRev')

    binRev = synthBinRev(rirHigh,psx.brir.binNoise,psx.brir.binRevBlockSizeRir,psx.brir.binRevBlockSizeNoise);

    %Shift in time according to peak of direct sound
    %Energy is smeared over time anyway, but a little time adjustment seams
    %reasonable
    [~,id] = max(abs(binRev));
    id = round(mean(id));

    if id < toaDirect
        binRev = [zeros(toaDirect-id,2);binRev];
        binRev = binRev(1:length(rir),:);
        binRev = binRev .* supdeq_win(length(binRev),[0 128]);
    elseif id > toaDirect
        binRev = binRev(id-toaDirect+1:end,:); 
        %Add zeros to have same length again
        binRev = [binRev;zeros(id-toaDirect,2)];
        binRev = binRev .* supdeq_win(length(binRev),[32 0]);
    end

    %Write in psx struct
    psx.brir.binRev = binRev;
    
else %If available, simply take from struct
    binRev = psx.brir.binRev;
end

%Filter if required
if hpBinRev
   [b,a] = butter(7,200/psx.fs*2,'high'); 
   binRev(:,1) = filter(b,a,binRev(:,1));
   binRev(:,2) = filter(b,a,binRev(:,2));
end

%% Find ID of BRIR where listener faces the source

if strcmp(sgSystem,'local')
    sgrad = sg * pi / 180;
    [xsg,ysg,zsg] = sph2cart(sgrad(:,1),sgrad(:,2),1);

    azelDirectRad = azelDirect * pi / 180;
    [xDirect,yDirect,zDirect] = sph2cart(azelDirectRad(:,1),azelDirectRad(:,2),1);

    %Calculate Euclidean distance
    euDistance = sqrt(sum(([xsg,ysg,zsg] - [xDirect,yDirect,zDirect]) .^ 2,2));
    %Get head orientation facing sound source and call frontalID
    [~,frontalID] = min(euDistance);
    %Write ID in struct
    psx.brir.frontalID = frontalID;
end

if strcmp(sgSystem,'global')
    sgrad = sg * pi / 180;
    [xsg,ysg,zsg] = sph2cart(sgrad(:,1),sgrad(:,2),1);

    azelDirectRad = [0,0]; %Just set to zero here instead of rotation. See below
    [xDirect,yDirect,zDirect] = sph2cart(azelDirectRad(:,1),azelDirectRad(:,2),1);

    %Calculate Euclidean distance
    euDistance = sqrt(sum(([xsg,ysg,zsg] - [xDirect,yDirect,zDirect]) .^ 2,2));
    %Get head orientation facing sound source and call frontalID
    [~,frontalID] = min(euDistance);
    %Write ID in struct
    psx.brir.frontalID = frontalID;
end

%% Get directivity for direct sound 
% This will always be used for compensation, as it is already part of the measurement

if applyDirectivity
    directivityMeas = AKisht(sourceDirectivity.SH.coeff, sourceDirectivity.SH.doFFT, [srcAzElDirect(1) 90-srcAzElDirect(2)], sourceDirectivity.SH.SHTmode, sourceDirectivity.SH.isEven, sourceDirectivity.SH.compact, sourceDirectivity.SH.SHmode);
    directivityMeas = real(directivityMeas);
    directivityMeas = AKinterpolateSpectrum(directivityMeas, sourceDirectivity.SH.f, directivityFilterLength, {'nearest' 'linear' 'nearest'}, sourceDirectivity.SH.fs);
    
    %Find nan if there is nan and replace with nearest neighbour
    m = flipud(directivityMeas);
    t = ~isnan(m);
    ii = cumsum(t);
    ii(ii == 0) = 1;
    ii = bsxfun(@plus,[0,ii(end,1:end-1)],ii);
    m1 = m(t);
    directivityMeas = flipud(m1(ii));
    clear m t ii m1;
    
    %Get directivity for extrapolated position too
    if strcmp(synthMode,'extrap')
        
        directivityEx = AKisht(sourceDirectivity.SH.coeff, sourceDirectivity.SH.doFFT, [srcAzElDirect_ex(1) 90-srcAzElDirect_ex(2)], sourceDirectivity.SH.SHTmode, sourceDirectivity.SH.isEven, sourceDirectivity.SH.compact, sourceDirectivity.SH.SHmode);
        directivityEx = real(directivityEx);
        directivityEx = AKinterpolateSpectrum(directivityEx, sourceDirectivity.SH.f, directivityFilterLength, {'nearest' 'linear' 'nearest'}, sourceDirectivity.SH.fs);

        %Find nan if there is nan and replace with nearest neighbour
        m = flipud(directivityEx);
        t = ~isnan(m);
        ii = cumsum(t);
        ii(ii == 0) = 1;
        ii = bsxfun(@plus,[0,ii(end,1:end-1)],ii);
        m1 = m(t);
        directivityEx = flipud(m1(ii));
        clear m t ii m1;
        
    end
end

%% Synthesize early part of BRIR

%Init variables
lengthHRIR = (length(hrtfs.f)*2-2)/hrtfs.FFToversize;

%Specular part of RIR
rirEarlySpec = rirHigh(1:mts2) .* wf;
%Diffuse part of RIR
rirEarlyDiff = binRev(1:mts2,:);

%Design inverse wf with a little reduced attenuation (no zero)
wfInv = sqrt(1-wf);
wfInvEnd = wfInv(directSE(2):end);
wfInvEnd(wfInvEnd < revAmplitude) = revAmplitude;
wfInv = [wfInv(1:directSE(2)-1);wfInvEnd];
rirEarlyDiff(:,1) = rirEarlyDiff(:,1).*wfInv;
rirEarlyDiff(:,2) = rirEarlyDiff(:,2).*wfInv;

%% SPAT mode
if strcmp(synthMode,'spat') %Mode to synthesize BRIR for spatialized measurement position
    
    disp('Synthesizing BRIRs...');
    
    %Define brirSpec as 3D matrix N x C x H, with N = number of samples, 
    %C = number of channels, H = number of head orientations
    brirSpec = zeros(mts2+lengthHRIR-1 + toaDirect,2,nHeadOr);

    %Rotate soundfield for global sampling grid system
    %Will results in (0,0) here, but just for the sake of it...
    if strcmp(sgSystem,'global')
        %[azelDirect_synth(1),azelDirect_synth(2)] = AKroomSimulationRotation(azelDirect(1),azelDirect(2),azelDirect(1),azelDirect(2));
        azelDirect_synth(1) = 0; azelDirect_synth(2) = 0;
        
        %Mirror sampling grid used for synthesis for global / relative 
        sg(:,1) = mod(-sg(:,1),360);
        sg(:,2) = sg(:,2)*-1;
        
    end
    if strcmp(sgSystem,'local') %No change
        azelDirect_synth = azelDirect;
    end
    
    %Synthesize direct sound first for each head orientation
    fprintf('|')
    for nho = 1:nHeadOr
        %Rotate again according to head orientations / sampling grid
        [azR,elR] = AKroomSimulationRotation(azelDirect_synth(1),azelDirect_synth(2),sg(nho,1),sg(nho,2));
        azelDirectSHrot = [azR,90-elR]; %Transform rotated elevation to SH coordinate format
        %Get HRIR
        [hrir(:,1),hrir(:,2)] = supdeq_getArbHRIR(hrtfs,azelDirectSHrot,'DEG',2,'ak');
        %Get IR part for filtering with soure directivity and hrir
        rirPart = rirEarlySpec(1:directSE(2));
        %Get rms of direct sound of rirEarlySpec/rirHigh
        eDirectEarly = rms(rirEarlySpec(directSE(1):directSE(2)));
        rirPartDirect = rirEarlySpec(directSE(1):directSE(2));
        
        if applyDirectivity
            %Apply source directivity
            %For direct sound of 'spat', this will simply be a dirac,
            %as the source directivity is already in the measurements, but this
            %way the procedure is the same for spat or extrap

            %df - directivity filter
            %df = directivityMeas ./ directivityMeas; --> dirac
            df = zeros(directivityFilterLength/2,1); %length/2 because of minimum phase filter to be used later
            df(1) = 1;
            %Filter direct sound with source directivity filter
            rirPart = convFFT(rirPart,df);
        end
        
        %Filter with HRIR
        dsBin = convFFT(rirPart,hrir);
        %Write in brirSpec array
        brirSpec(1:length(dsBin),:,nho) = brirSpec(1:length(dsBin),:,nho) + dsBin;
    end
    
    %Continue with early reflections
    %Rotate soundfield for global sampling grid system
    if strcmp(sgSystem,'global')
        [azelRef_synth(:,1),azelRef_synth(:,2)] = AKroomSimulationRotation(azelRef(:,1),azelRef(:,2),azelDirect(1),azelDirect(2));
    end
    if strcmp(sgSystem,'local') %No change
        azelRef_synth = azelRef;
    end
    %Synthesize early reflections for each head orientation
    for kk = 1:nRef
        fprintf('|')
        for nho = 1:nHeadOr
            %Rotate again according to head orientations / sampling grid
            [azR,elR] = AKroomSimulationRotation(azelRef_synth(kk,1),azelRef_synth(kk,2),sg(nho,1),sg(nho,2));
            azelRefSHrot = [azR,90-elR]; %Transform rotated elevation to SH coordinate format
            [hrir(:,1),hrir(:,2)] = supdeq_getArbHRIR(hrtfs,azelRefSHrot,'DEG',2,'ak');
            refBin = convFFT(rirEarlySpec(refSE(kk,1):refSE(kk,2)),hrir);
            brirSpec(refSE(kk,1):refSE(kk,2)+lengthHRIR-1,:,nho) = brirSpec(refSE(kk,1):refSE(kk,2)+lengthHRIR-1,:,nho) + refBin;
        end
    end  

    %Adjust TOA of frontal BRIR to directTOA
    toaBrirSpec = floor(mean(AKonsetDetect(squeeze(brirSpec(:,:,frontalID)),10,-20))); 
    brirSpec = circshift(brirSpec,toaDirect-toaBrirSpec);
    
    %Adjust brirSpec in level based on rms or magnitude of direct sound for head
    %orientation towards the source. This factor is required for all other
    %extrapolated BRIRs too.
    %Has to be done after circshift!
    if strcmp(drrMatch,'rms')
        rmsBrirSpecFrontalID = mean(rms(brirSpec(directSE(1):directSE(2),:,frontalID)));
        % Get compensation value based on eDirect of RIR parameterization
        brirSpecGain = eDirectEarly/rmsBrirSpecFrontalID;
    end
    if strcmp(drrMatch,'mag')
        %Do frequency dependent in reasonable frequency range based on
        %magnitude
        brirPartDirect = brirSpec(directSE(1):directSE(2),:,frontalID);
        NFFT = 2048;
        fVec = linspace(0,fs/2,NFFT/2+1);
        [~,fbin200] = min(abs(fVec-200));
        [~,fbin15k] = min(abs(fVec-15000));
        magRir = abs(fft(rirPartDirect,NFFT)); magRir = magRir(1:end/2+1);
        magBrir = abs(fft(brirPartDirect,NFFT)); magBrir = magBrir(1:end/2+1,:);
        M1 = mean(magRir(fbin200:fbin15k));
        B1 = mean(mean(magBrir(fbin200:fbin15k,:)));
        brirSpecGain = (M1/B1);
    end
    
    brirSpec = brirSpec*brirSpecGain;
    %Write in psx struct as required for extrapolated values
    psx.brir.brirSpecGain = brirSpecGain;
end

%% EXTRAP mode
if strcmp(synthMode,'extrap') %Mode to synthesize BRIR for extrapolated position

    disp('Synthesizing BRIRs...');
    
    %Define brirSpec as 3D matrix N x C x H, with N = number of samples, 
    %C = number of channels, H = number of head orientations
    brirSpec = zeros(mts2+lengthHRIR-1 + toaDirect,2,nHeadOr);
    
    %Rotate soundfield for global sampling grid system
    %Will results in (0,0) here, but just for the sake of it...
    if strcmp(sgSystem,'global')
        %[azelDirect__ex_synth(1),azelDirect_ex_synth(2)] = AKroomSimulationRotation(azelDirect_ex(1),azelDirect_ex(2),azelDirect_ex(1),azelDirect_ex(2));
        azelDirect_ex_synth(1) = 0; azelDirect_ex_synth(2) = 0;
        
        %Mirror sampling grid used for synthesis for global / relative 
        sg(:,1) = mod(-sg(:,1),360);
        sg(:,2) = sg(:,2)*-1;
        
    end
    if strcmp(sgSystem,'local') %No change
        azelDirect_ex_synth = azelDirect_ex;
    end
    %Synthesize direct sound first for each head orientation
    fprintf('|')
    for nho = 1:nHeadOr
        [azR,elR] = AKroomSimulationRotation(azelDirect_ex_synth(1),azelDirect_ex_synth(2),sg(nho,1),sg(nho,2));
        azelDirectSHrot = [azR,90-elR]; %Transform rotated elevation to SH coordinate format
        %Get HRIR
        [hrir(:,1),hrir(:,2)] = supdeq_getArbHRIR(hrtfs,azelDirectSHrot,'DEG',2,'ak');
        %Get IR part for filtering with soure directivity and hrir
        rirPart = rirEarlySpec(1:directSE(2));
        
        %Shift direct sound according to extrapolation
        rirShift = toaDirect_ex-toaDirect;
        if rirShift < 0 %toaDirect_ex < toaDirect %Decrease in distance, just do cirshift as its only zeros at the beginning
            rirPart = circshift(rirPart,rirShift);
        else %Increase in distance, add zeros
            rirPart = [zeros(rirShift,1);rirPart];
        end
        
        %Apply amplitude factor
        rirPart = rirPart*ampFactorDirect_ex;
        
        if applyDirectivity
            %Apply source directivity
            %For direct sound of 'spat', this will simply be a dirac,
            %as the source directivity is already in the measurements, but this
            %way the procedure is the same for spat or extrap

            %df - directivity filter
            df = directivityEx ./ directivityMeas;
            df(1) = 1; %Can be nan
            %Generate minimum phase filter based on magnitude response
            df = AKsingle2bothSidedSpectrum(df,true);
            df = real(ifft(df));
            df = circshift(df,directivityFilterLength/2-1);
            df = df .* hann(length(df));
            [~,df] = rceps(df);
            df = df(1:end/2,:);
            df = df.* supdeq_win(length(df),[0 length(df)/2]);
         
            %Filter direct sound with source directivity filter
            rirPart = convFFT(rirPart,df);
        end

        %Filter with HRIR
        dsBin = convFFT(rirPart,hrir);
        %Write in brirSpec array
        brirSpec(1:length(dsBin),:,nho) = brirSpec(1:length(dsBin),:,nho) + dsBin;
    end

    %Continue with early reflections
    %Rotate soundfield for global sampling grid system
    if strcmp(sgSystem,'global')
        [azelRef_ex_synth(:,1),azelRef_ex_synth(:,2)] = AKroomSimulationRotation(azelRef_ex(:,1),azelRef_ex(:,2),azelDirect(1),azelDirect(2));
    end
    if strcmp(sgSystem,'local') %No change
        azelRef_ex_synth = azelRef_ex;
    end
    %Synthesize early reflections for each head orientation
    for kk = 1:nRef
        fprintf('|')
        for nho = 1:nHeadOr
            [azR,elR] = AKroomSimulationRotation(azelRef_ex_synth(kk,1),azelRef_ex_synth(kk,2),sg(nho,1),sg(nho,2));
            azelRefSHrot = [azR,90-elR]; %Transform rotated elevation to SH coordinate format
            [hrir(:,1),hrir(:,2)] = supdeq_getArbHRIR(hrtfs,azelRefSHrot,'DEG',2,'ak');
            %Convolve with rir part of original RIR and shift later to new
            %position according to extrapolated values
            refBin = convFFT(rirEarlySpec(refSE(kk,1):refSE(kk,2)),hrir);
            %Apply corresponding amplitude factor
            refBin = refBin * refAmpFactor_ex(kk);
            %Place at new position in brirSpec array
            brirSpec(refSE_ex(kk,1):refSE_ex(kk,2)+lengthHRIR-1,:,nho) = brirSpec(refSE_ex(kk,1):refSE_ex(kk,2)+lengthHRIR-1,:,nho) + refBin;
        end
    end  

    %Adjust TOA of frontal BRIR to directTOA
    toaBrirSpec = floor(mean(AKonsetDetect(squeeze(brirSpec(:,:,frontalID)),10,-20))); 
    brirSpec = circshift(brirSpec,toaDirect_ex-toaBrirSpec);
    
    %Apply global gain value based on first spatialization
    brirSpec = brirSpec * brirSpecGain;
end

%% Combine specular and diffuse early part

%Shift rirEarlyDiff same as rirPart / direct sound in extrap mode
if strcmp(synthMode,'extrap')
   
    if rirShift < 0 %toaDirect_ex < toaDirect %Decrease in distance, just do cirshift as its only zeros at the beginning
        rirEarlyDiff = circshift(rirEarlyDiff,rirShift);
    else %Increase in distance, add zeros
        rirEarlyDiff = [zeros(rirShift,2);rirEarlyDiff];
    end
    
end

brirEarly = zeros(length(rirEarlyDiff),2,nHeadOr);
for nho = 1:nHeadOr
    brirEarly(1:length(rirEarlyDiff),:,nho) = brirSpec(1:length(rirEarlyDiff),:,nho) + rirEarlyDiff;
end

%% Combine parts to entire BRIR

%Shift bin same as rirPart / direct sound in extrap mode
if strcmp(synthMode,'extrap')

    if rirShift < 0 %toaDirect_ex < toaDirect %Decrease in distance, just do cirshift as its only zeros at the beginning
        binRev = circshift(binRev,rirShift);
    else %Increase in distance, add zeros
        binRev = [zeros(rirShift,2);binRev];
    end
    
end

%Combine early and late part
synthBRIR = zeros(size(binRev,1),size(binRev,2),nHeadOr);
for nho = 1:nHeadOr
    synthBRIR(:,:,nho) = binRev;
    synthBRIR(1:size(brirEarly,1),:,nho) = brirEarly(:,:,nho);
end

%Add low frequency component again.
if applyXover
    
    %Shift rirLow and adjust level in extrap mode
    if strcmp(synthMode,'extrap')
        if rirShift < 0 %toaDirect_ex < toaDirect %Decrease in distance, just do cirshift as its only zeros at the beginning
            rirLow = circshift(rirLow,rirShift);
        else %Increase in distance, add zeros
            rirLow = [zeros(rirShift,1);rirLow];
        end
        rirLow = rirLow*ampFactorDirect_ex;
    end
    
    rirLow = rirLow * brirSpecGain;
    
    for nho = 1:nHeadOr
        synthBRIR(:,1,nho)=synthBRIR(:,1,nho)+rirLow;
        synthBRIR(:,2,nho)=synthBRIR(:,2,nho)+rirLow;
    end
end

%% Apply filtering 

%Interaural coherence filter
if applyIC
   for nho = 1:nHeadOr
       synthBRIR(:,:,nho) = AKbinauralCoherence(synthBRIR(:,:,nho),fs);
   end
end

%Measurement filter
if applyMeasFilt
    synthBRIRFilt = zeros(size(synthBRIR,1)+length(measFiltMin)-1,size(synthBRIR,2),size(synthBRIR,3));
    for nho = 1:nHeadOr
        synthBRIRFilt(:,:,nho) = convFFT(measFiltMin,synthBRIR(:,:,nho));
    end
    synthBRIR = synthBRIRFilt;
end

%Nth order Butterworth lowcut/highpass
if butterFc > 0
    [bLowcut,aLowcut] = butter(butterN,butterFc/fs*2,'high'); 
    for nho = 1:nHeadOr
        synthBRIR(:,:,nho) = filter(bLowcut,aLowcut,synthBRIR(:,:,nho));
    end
end

%% Save in struct

if strcmp(synthMode,'spat')
    psx.brir.synthBRIR_spat = synthBRIR;
end

if strcmp(synthMode,'extrap')
    psx.brir.synthBRIR_extrap = synthBRIR;
end

fprintf('\nDone...\n');

end