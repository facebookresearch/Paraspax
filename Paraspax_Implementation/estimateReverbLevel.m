% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function psx = estimateReverbLevel(psx,relToDS)
%
% Functions to estimate the reverberation level in a monaural IR 
% based on different estimation approaches. In general, the reverberation
% level at the mixing time according to Abel is determined and from this
% point different slopes are applied. Another approach simply determines
% the residual energy after substracting the specular reflections with the
% weighting function generated in the function "detectReflections".
%
% Output:
% psx                   - psx struct with reverberation level (based on
%                         EDC) and different reverberation level onset 
%                         points in field psx.rev.
%
% Input:        
% psx                   - psx struct with required fields
% relToDS               - Boolean to decide whether the levels should be
%                         determined as abolute levels (relToDS = false) or 
%                         relative to the maximum absolute level of the 
%                         direct sound (relToDS = true)
%                         Default: true (in most cases better as values are
%                         used to design weighting function)
%
% Dependencies: -
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function psx = estimateReverbLevel(psx,relToDS)

if nargin < 2
    relToDS = true;
end

%Get various variables for further processing
fs  = psx.fs;
rir = psx.rir;
tDirect = psx.par.mon.rirOnset_20dB; %TOA of direct sound in samples
directWin = round(psx.refDet.directWin*fs); %Direct sound window in samples
refWin = round(psx.refDet.refWin*fs); %Direct sound window in samples
eDirect = psx.refDet.eDirect;
mts = round(psx.par.mon.mtAbel * 10^-3 * fs); %Mixing time according to abel in samples
wf  = psx.refDet.wf; %Weighting function
estMethod = psx.rev.method; %Estimation method
dsl = max(abs(rir(tDirect-directWin(1):tDirect+directWin(2))));%Get level of direct sound

%% Get info about first reflection

[tFirstRef,firstRefID] = min(psx.refDet.refListDet(:,1));
eFirstRef = psx.refDet.refListDet(firstRefID,2);

%% EDC method

if strcmp(estMethod,'EDC') || strcmp(estMethod,'EDC_RMS_MAX')
    
    %Transform EDC (energy curve) to level curve
    edc = psx.par.mon.EDC_BB_AK; %Broadband EDC estimated with AKtools
    edc = 10.^(edc/10);
    edc = sqrt(edc);
    
    %Shift edc (which is normalized in level to rir)
    directTaps = tDirect-directWin(1):tDirect+directWin(2);
    eEDCdirect = rms(edc(directTaps));
    edc = edc * (eDirect*eEDCdirect);

    %%Get level at mixing time
    mtWin = round(0.001*fs);
    if mod(mtWin,2)
        mtWin = mtWin+1;
    end
    ampAtMT = mean(edc(mts-mtWin/2+1:mts+mtWin/2));
    
    %Get different estimates
    %(1) - Straight line from start of RIR to mixing time with amplitude of
    %mixing time, after that level decrease according to edc
    straightLine = ones(length(rir),1)*ampAtMT;
    straightLine(mts+1:end) = edc(mts+1:end);
    straightLine = straightLine.';
    %Reverberation amplitude at TOA of direct sound and first reflection (Will
    %be the same obviously...)
    sLine_revAmp_direct = straightLine(tDirect);
    sLine_revAmp_firstRef = straightLine(tFirstRef);

    %(2) - Logarithmic level increase (linear in amplitude) from first
    %reflexion to mixing time. Start 60dB below level at mixing time.
    tVec = (tFirstRef:mts);
    startVal = eFirstRef * 10^-3;%-60dB below energy of first reflection
    logLine = linspace(startVal,ampAtMT,length(tVec)); %Logarithmic in dB, Linear in amplitude
    logLine = [ones(1,tFirstRef-1)*startVal,logLine,edc(mts+1:end)'];
    %Reverberation amplitude at TOA of direct sound and first reflection
    logLine_revAmp_direct = logLine(tDirect);
    logLine_revAmp_firstRef = logLine(tFirstRef);

    %(3) - Linear level increase (exponential in amplitude) from first
    %reflexion to mixing time. Start 60dB below level at mixing time.
    linLine = exp(linspace(log(startVal),log(ampAtMT),length(tVec)));
    linLine = [ones(1,tFirstRef-1)*startVal,linLine,edc(mts+1:end)'];
    %Reverberation amplitude at TOA of direct sound and first reflection
    linLine_revAmp_direct = linLine(tDirect);
    linLine_revAmp_firstRef = linLine(tFirstRef);

    %(4) - Linear interpolation of edc level curve from mixing time to start
    x = 2*mts:3*mts;
    interpLine = edc(x(1):x(end))';
    interpLine = 20*log10(interpLine);
    p = polyfit(x,interpLine,1);
    f = polyval(p,x);
    interpLine = interp1(x(1):x(end),f,1:length(edc),'linear','extrap');
    interpLine = 10.^(interpLine/20);
    %Reverberation amplitude at TOA of direct sound and first reflection
    interpLine_revAmp_direct = interpLine(tDirect);
    interpLine_revAmp_firstRef = interpLine(tFirstRef);
    interpLine_revAmp_mt = interpLine(mts);
    
    %Adjust relative to direct sound level if desired
    if relToDS
        edc = edc/dsl;
        ampAtMT = ampAtMT/dsl;
        straightLine = straightLine/dsl;
        sLine_revAmp_direct = sLine_revAmp_direct/dsl;
        sLine_revAmp_firstRef = sLine_revAmp_firstRef/dsl;
        logLine = logLine/dsl;
        logLine_revAmp_direct = logLine_revAmp_direct/dsl;
        logLine_revAmp_firstRef = logLine_revAmp_firstRef/dsl;
        linLine = linLine/dsl;
        linLine_revAmp_direct = linLine_revAmp_direct/dsl;
        linLine_revAmp_firstRef = linLine_revAmp_firstRef/dsl;
        interpLine = interpLine/dsl;
        interpLine_revAmp_direct = interpLine_revAmp_direct/dsl;
        interpLine_revAmp_firstRef = interpLine_revAmp_firstRef/dsl;
        interpLine_revAmp_mt = interpLine_revAmp_mt/dsl;
    end
    
    %Write results in struct
    psx.rev.tFirstRef = tFirstRef; %TOA of first reflection
    psx.rev.relToDS = relToDS; %Info if relative or absolute levels
    psx.rev.edcMethod.edcLevel = edc; %EDC as adjusted level curve of reverberation
    psx.rev.edcMethod.ampAtMT = ampAtMT; %Amplitude at mixing time
    psx.rev.edcMethod.straightLine = straightLine; %Reverberation level with straight line from mixing time
    psx.rev.edcMethod.straightLine_revAmp_direct = sLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.edcMethod.straightLine_revAmp_firstRef = sLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.edcMethod.logLine = logLine; %Reverberation level with log line from mixing time
    psx.rev.edcMethod.logLine_revAmp_direct = logLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.edcMethod.logLine_revAmp_firstRef = logLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.edcMethod.linLine = linLine; %Reverberation level with linear line from mixing time
    psx.rev.edcMethod.linLine_revAmp_direct = linLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.edcMethod.linLine_revAmp_firstRef = linLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.edcMethod.interpLine = interpLine; %Reverberation level with interpolated line from mixing time
    psx.rev.edcMethod.interpLine_revAmp_direct = interpLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.edcMethod.interpLine_revAmp_firstRef = interpLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.edcMethod.interpLine_revAmp_mt = interpLine_revAmp_mt; %Amplitude at MT

end

%% RMS method

if strcmp(estMethod,'RMS') || strcmp(estMethod,'EDC_RMS_MAX')
    
    rmsCurve = slidingRMS(rir,round(0.001*fs));
    
    %%Get level at mixing time
    mtWin = round(0.001*fs);
    if mod(mtWin,2)
        mtWin = mtWin+1;
    end
    ampAtMT = mean(rmsCurve(mts-mtWin/2+1:mts+mtWin/2));
    
    %Get different estimates
    %Get polyfit of RMS curve first 
    %(1) - Linear interpolation of rms level curve from mixing time to start
    x = 2*mts:3*mts;
    interpLine = rmsCurve(x(1):x(end))';
    interpLine = 20*log10(interpLine);
    p = polyfit(x,interpLine,1);
    f = polyval(p,x);
    interpLine = interp1(x(1):x(end),f,1:length(rmsCurve),'linear','extrap');
    interpLine = 10.^(interpLine/20);
    %Reverberation amplitude at TOA of direct sound and first reflection
    interpLine_revAmp_direct = interpLine(tDirect);
    interpLine_revAmp_firstRef = interpLine(tFirstRef);
    interpLine_revAmp_mt = interpLine(mts);
    
    %(2) - Straight line from start of RIR to mixing time with amplitude of
    %mixing time, after that level decrease according to edc
    straightLine = ones(length(rir),1)*ampAtMT;
    straightLine(mts+1:end) = interpLine(mts+1:end);
    %Reverberation amplitude at TOA of direct sound and first reflection (Will
    %be the same obviously...)
    sLine_revAmp_direct = straightLine(tDirect);
    sLine_revAmp_firstRef = straightLine(tFirstRef);

    %(2) - Logarithmic level increase (linear in amplitude) from first
    %reflexion to mixing time. Start 60dB below level at mixing time.
    tVec = (tFirstRef:mts);
    startVal = eFirstRef * 10^-3;%-60dB below energy of first reflection
    logLine = linspace(startVal,ampAtMT,length(tVec)); %Logarithmic in dB, Linear in amplitude
    logLine = [ones(1,tFirstRef-1)*startVal,logLine,interpLine(mts+1:end)];
    %Reverberation amplitude at TOA of direct sound and first reflection
    logLine_revAmp_direct = logLine(tDirect);
    logLine_revAmp_firstRef = logLine(tFirstRef);

    %(3) - Linear level increase (exponential in amplitude) from first
    %reflexion to mixing time. Start 60dB below level at mixing time.
    linLine = exp(linspace(log(startVal),log(ampAtMT),length(tVec)));
    linLine = [ones(1,tFirstRef-1)*startVal,linLine,interpLine(mts+1:end)];
    %Reverberation amplitude at TOA of direct sound and first reflection
    linLine_revAmp_direct = linLine(tDirect);
    linLine_revAmp_firstRef = linLine(tFirstRef);
    
    %Adjust relative to direct sound level if desired
    if relToDS
        rmsCurve = rmsCurve/dsl;
        ampAtMT = ampAtMT/dsl;
        straightLine = straightLine/dsl;
        sLine_revAmp_direct = sLine_revAmp_direct/dsl;
        sLine_revAmp_firstRef = sLine_revAmp_firstRef/dsl;
        logLine = logLine/dsl;
        logLine_revAmp_direct = logLine_revAmp_direct/dsl;
        logLine_revAmp_firstRef = logLine_revAmp_firstRef/dsl;
        linLine = linLine/dsl;
        linLine_revAmp_direct = linLine_revAmp_direct/dsl;
        linLine_revAmp_firstRef = linLine_revAmp_firstRef/dsl;
        interpLine = interpLine/dsl;
        interpLine_revAmp_direct = interpLine_revAmp_direct/dsl;
        interpLine_revAmp_firstRef = interpLine_revAmp_firstRef/dsl;
        interpLine_revAmp_mt = interpLine_revAmp_mt/dsl;
    end
    
    %Write results in struct
    psx.rev.tFirstRef = tFirstRef; %TOA of first reflection
    psx.rev.relToDS = relToDS; %Info if relative or absolute levels
    psx.rev.rmsMethod.rmsCurve = rmsCurve; %RMS curve obtained with sliding window rms
    psx.rev.rmsMethod.ampAtMT = ampAtMT; %Amplitude at mixing time
    psx.rev.rmsMethod.straightLine = straightLine; %Reverberation level with straight line from mixing time
    psx.rev.rmsMethod.straightLine_revAmp_direct = sLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.rmsMethod.straightLine_revAmp_firstRef = sLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.rmsMethod.logLine = logLine; %Reverberation level with log line from mixing time
    psx.rev.rmsMethod.logLine_revAmp_direct = logLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.rmsMethod.logLine_revAmp_firstRef = logLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.rmsMethod.linLine = linLine; %Reverberation level with linear line from mixing time
    psx.rev.rmsMethod.linLine_revAmp_direct = linLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.rmsMethod.linLine_revAmp_firstRef = linLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.rmsMethod.interpLine = interpLine; %Reverberation level with interpolated line from mixing time
    psx.rev.rmsMethod.interpLine_revAmp_direct = interpLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.rmsMethod.interpLine_revAmp_firstRef = interpLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.rmsMethod.interpLine_revAmp_mt = interpLine_revAmp_mt; %Amplitude at MT
    
end

%% MAX method

if strcmp(estMethod,'MAX') || strcmp(estMethod,'EDC_RMS_MAX')
    
    maxCurve = slidingMax(rir,round(0.001*fs));
    
    %%Get level at mixing time
    mtWin = round(0.001*fs);
    if mod(mtWin,2)
        mtWin = mtWin+1;
    end
    ampAtMT = mean(maxCurve(mts-mtWin/2+1:mts+mtWin/2));
    
    %Get different estimates
    %Get polyfit of sliding max curve first 
    %(1) - Linear interpolation of sliding max curve from mixing time to start
    x = 2*mts:3*mts;
    interpLine = maxCurve(x(1):x(end))';
    interpLine = 20*log10(interpLine);
    p = polyfit(x,interpLine,1);
    f = polyval(p,x);
    interpLine = interp1(x(1):x(end),f,1:length(maxCurve),'linear','extrap');
    interpLine = 10.^(interpLine/20);
    %Reverberation amplitude at TOA of direct sound and first reflection
    interpLine_revAmp_direct = interpLine(tDirect);
    interpLine_revAmp_firstRef = interpLine(tFirstRef);
    interpLine_revAmp_mt = interpLine(mts);
    
    %(2) - Straight line from start of RIR to mixing time with amplitude of
    %mixing time, after that level decrease according to edc
    straightLine = ones(length(rir),1)*ampAtMT;
    straightLine(mts+1:end) = interpLine(mts+1:end);
    %Reverberation amplitude at TOA of direct sound and first reflection (Will
    %be the same obviously...)
    sLine_revAmp_direct = straightLine(tDirect);
    sLine_revAmp_firstRef = straightLine(tFirstRef);

    %(2) - Logarithmic level increase (linear in amplitude) from first
    %reflexion to mixing time. Start 60dB below level at mixing time.
    tVec = (tFirstRef:mts);
    startVal = eFirstRef * 10^-3;%-60dB below energy of first reflection
    logLine = linspace(startVal,ampAtMT,length(tVec)); %Logarithmic in dB, Linear in amplitude
    logLine = [ones(1,tFirstRef-1)*startVal,logLine,interpLine(mts+1:end)];
    %Reverberation amplitude at TOA of direct sound and first reflection
    logLine_revAmp_direct = logLine(tDirect);
    logLine_revAmp_firstRef = logLine(tFirstRef);

    %(3) - Linear level increase (exponential in amplitude) from first
    %reflexion to mixing time. Start 60dB below level at mixing time.
    linLine = exp(linspace(log(startVal),log(ampAtMT),length(tVec)));
    linLine = [ones(1,tFirstRef-1)*startVal,linLine,interpLine(mts+1:end)];
    %Reverberation amplitude at TOA of direct sound and first reflection
    linLine_revAmp_direct = linLine(tDirect);
    linLine_revAmp_firstRef = linLine(tFirstRef);
    
    %Adjust relative to direct sound level if desired
    if relToDS
        maxCurve = maxCurve/dsl;
        ampAtMT = ampAtMT/dsl;
        straightLine = straightLine/dsl;
        sLine_revAmp_direct = sLine_revAmp_direct/dsl;
        sLine_revAmp_firstRef = sLine_revAmp_firstRef/dsl;
        logLine = logLine/dsl;
        logLine_revAmp_direct = logLine_revAmp_direct/dsl;
        logLine_revAmp_firstRef = logLine_revAmp_firstRef/dsl;
        linLine = linLine/dsl;
        linLine_revAmp_direct = linLine_revAmp_direct/dsl;
        linLine_revAmp_firstRef = linLine_revAmp_firstRef/dsl;
        interpLine = interpLine/dsl;
        interpLine_revAmp_direct = interpLine_revAmp_direct/dsl;
        interpLine_revAmp_firstRef = interpLine_revAmp_firstRef/dsl;
        interpLine_revAmp_mt = interpLine_revAmp_mt/dsl;
    end
    
    %Write results in struct
    psx.rev.tFirstRef = tFirstRef; %TOA of first reflection
    psx.rev.relToDS = relToDS; %Info if relative or absolute levels
    psx.rev.maxMethod.maxCurve = maxCurve; %RMS curve obtained with sliding window rms
    psx.rev.maxMethod.ampAtMT = ampAtMT; %Amplitude at mixing time
    psx.rev.maxMethod.straightLine = straightLine; %Reverberation level with straight line from mixing time
    psx.rev.maxMethod.straightLine_revAmp_direct = sLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.maxMethod.straightLine_revAmp_firstRef = sLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.maxMethod.logLine = logLine; %Reverberation level with log line from mixing time
    psx.rev.maxMethod.logLine_revAmp_direct = logLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.maxMethod.logLine_revAmp_firstRef = logLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.maxMethod.linLine = linLine; %Reverberation level with linear line from mixing time
    psx.rev.maxMethod.linLine_revAmp_direct = linLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.maxMethod.linLine_revAmp_firstRef = linLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.maxMethod.interpLine = interpLine; %Reverberation level with interpolated line from mixing time
    psx.rev.maxMethod.interpLine_revAmp_direct = interpLine_revAmp_direct; %Amplitude at tDirect
    psx.rev.maxMethod.interpLine_revAmp_firstRef = interpLine_revAmp_firstRef; %Amplitude at tRef
    psx.rev.maxMethod.interpLine_revAmp_mt = interpLine_revAmp_mt; %Amplitude at MT
    
end
end