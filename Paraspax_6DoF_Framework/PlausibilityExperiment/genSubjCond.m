% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function genSubjCond()
%
% Function to randomize the conditions of the 6-DoF plausibility experiment 
% for the respective participant. For more details on the Paraspax real-time 
% framework and the plausibility study, see [1].
%
% Output:
% -                   
%
% Input:        
% -
%
% Dependencies: oscsend, global variable 'sh' and 'binsim' providing
% information about udp and cartesian/spherical grid, Natnetbib from
% OptiTrack, Paraspax real-time framework, PyBinSim, OptiTrack
%
% References:
% [1] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on One Monaural Room Impulse Response,” 
% J. Audio Eng. Soc., vol. 69, no. 7/8, pp. 557–575, 2021. ﻿
% https://doi.org/10.17743/jaes.2021.0009
%
% Code written 2019/2020 by JMA, Johannes M. Arend.


function genSubjCond()

global subj;
global sh;
global binsim;

subj.expDate = datetime('now');
subj.trial = 0; %Counter for trials

%Calculate real/virtual sequence (Not exactly 50/50 to reduce
%predictability, +- 10% of prop). EDIT:CHANGED
nTrials = sh.nTestSignals * sh.nSources;
zeroNum = zeros(nTrials,1);
oneNum = ones(nTrials,1);
prop = nTrials / 2;
%perc = ceil(prop * 0.10);
%prop = prop + (randi([-perc,perc],1)); %Changed to 50/50 as not necessary
%to slighty vary numbe of real/virtual sources. Participants are not able
%to count anyway after the first insecure answer...
sequ = [zeroNum(1:prop);oneNum(1:abs(nTrials-prop))];
sequ = sequ(randperm(numel(sequ)));

%Write random test signal - source vector
testSigSequ = 1:sh.nTestSignals;
testSigSequ = repmat(testSigSequ,1,sh.nSources); testSigSequ = testSigSequ';
sourceSequ  = repmat(binsim.sourceIDs,sh.nTestSignals,1); sourceSequ = sourceSequ(:);
%Randomize with same random permuation
randPer = randperm(nTrials);
testSigSequ = testSigSequ(randPer);
sourceSequ = sourceSequ(randPer);

subj.nTrials = nTrials;
if isfield(sh,'nPositions')
    subj.nTrialsPerPosition = subj.nTrials / sh.nPositions;
    if floor(subj.nTrialsPerPosition) ~= subj.nTrialsPerPosition
        error('Please provide appropriate number of trials and positions');
    end
end
subj.rvsequ = sequ; %Real/Virtual sequence; 0 - Real, 1 - Virtual
subj.testSigSequ = testSigSequ; %Randomized sequence of test signals
subj.sourceSequ = sourceSequ; %Randomized sequence of source IDs
subj.prop = [prop, abs(nTrials-prop)];
subj.run = 1; %Save that already started
subj.resp = nan(nTrials,1); %For raw response if choosen real / virtual (left / right on controller)
subj.respCorrect = nan(nTrials,1); %For correct responses
subj.expComp = 0; %Variable to say if experiment was completed

%Prepare logfile
currentTime = clock;
timeStr = [num2str(currentTime(4)),'-',num2str(currentTime(5))];
subj.logtxt = fopen(fullfile('subjects',[subj.ID,'_',date,'_',timeStr,'_log.txt']),'a');
if subj.logtxt == -1
   error('Cannot open log file!'); 
end

subj.fileName = ['subjects/',subj.ID,'_',date,'_',timeStr];
save(subj.fileName,'subj');