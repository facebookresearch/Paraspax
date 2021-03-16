% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function rir_post = postProRIR(rir,fs,fsTarget,normValue,sourceDistance,doa)
%
% This function applies basic post processing to raw measurement data.
% Input can be single channel (monaural RIR) or two channel (BRIR).
%
% Output:
% rir_post              - Processed rir
% doa_post              - Processed (time-shifted) DOA vector if passed
%
% Input:        
% rir                   - Raw RIR
% fs                    - Sample rate of raw RIR
% fsTarget              - Target sample rate of post processed rir if
%                         up/downsampling should be applied. Leave empty if
%                         sample rate should not be changed
%                         Default: no sample rate conversion
% normValue             - Apply peak normalization if required. Define
%                         peak value (linear) if normalization should be
%                         applied. If set to 0, no normalization will be
%                         applied
%                         Default: 0 (no normalization)
% sourceDistance        - Source distance in m
%                         If defined, TOA of RIR will be shifted according 
%                         to given source distance. Set to 0 if no shift of
%                         RIR should be applied
%                         Default: 0 (no shift)
% doa                   - If DOA vector from SDM analysis is passed, RIR is
%                         only normalized and shifted according to source
%                         distance (better for extrapolation) and DOA
%                         vector is shifted in the same way to keep
%                         relation between RIR and DOA vector 
%                         If empty, normal pre-processing is gonne be
%                         applied
%
%
% Dependencies: SUpDEq toolbox, AKtools, ITA toolbox
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function [rir_post,doa_post] = postProRIR(rir,fs,fsTarget,normValue,sourceDistance,doa)

if size(rir,2) > size(rir,1)
    rir = rir.';
end

if nargin < 2
    error('Please specify sample rate of RIR');
end

if nargin < 3 || isempty(fsTarget)
    fsTarget = fs; % no resampling;
end

if nargin < 4 || isempty(normValue)
    normValue = 0; % no peak normalization;
end

if nargin < 5 || isempty(sourceDistance)
    sourceDistance = 0; % no shift
end

if nargin < 6 || isempty(doa)
    doa = []; % no DOA vector related processing
end

%Check if two channel RIR (BRIR) was passed
brir = false;
if size(rir,2) == 2
    brir = true;
end

%% Post-processing

if isempty(doa)

    %Get onset 
    rirOnset = floor(AKonsetDetect(rir,10,-20,'rel')); %Set to -20dB in accordance with ITA toolbox
    %Choose smallest value if brir
    if brir
        rirOnset = min(rirOnset);
    end
    %Shift 5 ms
    rirCutHead = rirOnset - ceil(0.005*fs);
    if rirCutHead < 1
        rirCutHead = 1;
    end
    %Cut head
    rir_post = rir(rirCutHead:end,:);

    %Get transition to noise floor (broadband) with lundeby noise estimation in
    %samples
    rir_post_ITA = var2ITA(rir_post,fs); %Needs to be an ITA object
    rirCutTail = ita_roomacoustics(rir_post_ITA,'broadbandAnalysis','Intersection_Time_Lundeby');
    rirCutTail = ceil(rirCutTail.Intersection_Time_Lundeby.freqData * fs);
    %Choose smallest value if brir
    if brir
        rirCutTail = min(rirCutTail);
    end
    %Shift 5 ms
    rirCutTail = rirCutTail+0.005*fs;
    %Cut tail
    rir_post = rir_post(1:rirCutTail,:);

    %Apply 5 ms window to head and tail
    win = supdeq_win(length(rir_post),[floor(0.005*fs)-1,floor(0.005*fs)-1]);
    if brir
        rir_post(:,1) = rir_post(:,1).*win;
        rir_post(:,2) = rir_post(:,2).*win;
    else
        rir_post = rir_post.*win;
    end

    %Apply resampling if desired
    if fs ~= fsTarget
        gComDiv = gcd(fs, fsTarget);
        p = fsTarget / gComDiv;
        q = fs / gComDiv;
        rir_post = resample(rir_post,double(p),double(q));
    end

    %Apply peak normalization if desired
    if normValue ~= 0
        if brir
            rir_post = rir_post * normValue/max(max(abs(rir_post))); 
        else
            rir_post = rir_post * normValue/max(abs(rir_post)); 
        end
    end

    %Shift according to source distance if desired
    if sourceDistance ~= 0
        rirOnset = floor(min(AKonsetDetect(rir_post,10,-20,'rel'))); %Set to -20dB in accordance with ITA toolbox
        tDirect = round(sourceDistance*fsTarget/343);
        if rirOnset < tDirect
            if brir
                rir_post = [zeros(tDirect-rirOnset,2);rir_post];
            else
                rir_post = [zeros(tDirect-rirOnset,1);rir_post];
            end
        end
        
        if rirOnset > tDirect
            rir_post = rir_post(floor(rirOnset-tDirect-1):end,:);
        end

    end

end

% DOA processing if DOA vector passed
if ~isempty(doa)
    
    %Normalization if desired
    if normValue ~= 0
        if brir
            rir_post = rir * normValue/max(max(abs(rir))); 
        else
            rir_post = rir * normValue/max(abs(rir)); 
        end
    end

    %Shift according to source distance if desired
    if sourceDistance ~= 0
        rirOnset = floor(min(AKonsetDetect(rir_post,10,-20,'rel'))); %Set to -20dB in accordance with ITA toolbox
        tDirect = round(sourceDistance*fsTarget/343);
        if rirOnset < tDirect
            if brir
                rir_post = [zeros(tDirect-rirOnset,2);rir_post];
            else
                rir_post = [zeros(tDirect-rirOnset,1);rir_post];
            end
            doa_post = [zeros(tDirect-rirOnset,3);doa];
        end

        if rirOnset > tDirect
            rir_post = rir_post(floor(rirOnset-tDirect-1):end,:);
            doa_post = doa_post(floor(rirOnset-tDirect-1):end,:);
        end

    end
    
end