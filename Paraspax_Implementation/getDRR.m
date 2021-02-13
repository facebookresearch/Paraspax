%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function DRR = getDRR(rir,winLength,fs)
%
% Function to calculate direct-to-reverberation ratio (DRR) in dB
%
% Output:
% DRR                   - DRR value in dB
%
% Input:        
% rir                   - Room impulse response. If two channel IR (BRIR)
%                         is passed, the combined/binaural DRR is
%                         calculated
% winLength             - Window length (in ms) assigned to direct sound.
%                         Window will be +winLength after detected onset
%                         of direct sound (with safety margin before direct
%                         sound of 0.5 ms)
%                         Default: 1 ms
% fs                    - Sample rate of rir
%                         Default: 48 kHz
%
%
% Dependencies: AKtools
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.
%
% Copyright (c) Facebook, Inc. and its affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

function DRR = getDRR(rir,winLength,fs)

if size(rir,2) > size(rir,1)
    rir = rir.';
end

brir = false;
if size(rir,2) == 2
    brir = true;
end

if nargin < 2
    winLength = 1;
end

if nargin < 3 || isempty(fs)
    fs = 48000;
end

%% Calculate broadband DRR

if ~brir %Single-Channel RIR
    
    rirOnset = floor(AKonsetDetect(rir,10,-20,'rel')); %Set to -20dB in accordance with ITA toolbox and Zahorik2002

    cplus = round((winLength/1000)*fs); %Paraspax according to window later used for energy detection in direct sound
    cminus = round(0.0005*fs); %Safety margin of 0.5 ms

    rirDir = rir(rirOnset-cminus:rirOnset+cplus); 
    Edir = sum(abs(rirDir).^2);

    rirRev = rir(rirOnset+cplus+1:end);
    Erev = sum(abs(rirRev).^2);

    DRR = 10*log10(Edir / Erev);
    
else %Two-Channel RIR / BRIR
    
    rirOnsetL = floor(AKonsetDetect(rir(:,1),10,-20,'rel'));
    rirOnsetR = floor(AKonsetDetect(rir(:,2),10,-20,'rel'));
    
    cplus = round((winLength/1000)*fs);
    cminus = round(0.0005*fs);
    
    rirDirL = rir(rirOnsetL-cminus:rirOnsetL+cplus,1); 
    EdirL = sum(abs(rirDirL).^2);
    rirDirR = rir(rirOnsetR-cminus:rirOnsetR+cplus,2); 
    EdirR = sum(abs(rirDirR).^2);
    
    rirRevL = rir(rirOnsetL+cplus+1:end,1);
    ErevL = sum(abs(rirRevL).^2);
    rirRevR = rir(rirOnsetR+cplus+1:end,2);
    ErevR = sum(abs(rirRevR).^2);
    
    DRR = 10*log10( (EdirL + EdirR) / (ErevL + ErevR) );

end