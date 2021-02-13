%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function genSubjCond_Positions()
%
% Function to randomize the conditions of the 3-DoF plausibility experiment 
% for the respective participant. For more details on the Paraspax real-time 
% frameweork and the Plausibility study, see [1]. 'Positions' is for a
% 3-DoF plausibility study.
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

function genSubjCond_Positions()

global subj;
global sh;
global binsim;

subj.expDate = datetime('now');
subj.trial = 0; %Counter for trials
subj.currentListenerPosition = 1; %Counter for positions to be tested
subj.expCompletePerPosition = zeros(sh.nPositions,1);

%Calculate real/virtual sequence (Not exactly 50/50 to reduce
%predictability, +- 10% of prop). EDIT:CHANGED
nTrials = sh.nTestSignals * sh.nSources;
nTrialsPerPosition = sh.nTestSignals * sh.nSources / sh.nPositions;
if floor(nTrialsPerPosition) ~= nTrialsPerPosition
    error('Please provide appropriate number of trials and positions');
end
zeroNum = zeros(nTrialsPerPosition,1);
oneNum = ones(nTrialsPerPosition,1);
prop = nTrialsPerPosition / 2;
%perc = ceil(prop * 0.10);
%prop = prop + (randi([-perc,perc],1)); %Changed to 50/50 as not necessary
%to slighty vary numbe of real/virtual sources. Participants are not able
%to count anyway after the first insecure answer...
sequ = [zeroNum(1:prop);oneNum(1:abs(nTrialsPerPosition-prop))];
for kk = 1:sh.nPositions
    sequPerPosition(:,kk) = sequ(randperm(numel(sequ)));
end
sequPerPosition = sequPerPosition(:);

%Write random test signal - source vector
testSigSequ = 1:sh.nTestSignals;
testSigSequ = repmat(testSigSequ,1,sh.nSources); testSigSequ = testSigSequ';
sourceSequ  = repmat(binsim.sourceIDs,sh.nTestSignals,1); sourceSequ = sourceSequ(:);
%Randomize with same random permuation
randPer = randperm(nTrials);
testSigSequ = testSigSequ(randPer);
sourceSequ = sourceSequ(randPer);


subj.nTrials = nTrials;
subj.nTrialsPerPosition = nTrialsPerPosition; %Trials per listener position
subj.rvsequ = sequPerPosition; %Real/Virtual sequence per position; 0 - Real, 1 - Virtual
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
%Open later in position experiment
%subj.logtxt = fopen(fullfile('subjects',[subj.ID,'_',date,'_',timeStr,'_log.txt']),'a');
%if subj.logtxt == -1
%   error('Cannot open log file!'); 
%end
subj.logtxtName = [subj.ID,'_Position','_',date,'_',timeStr,'_log.txt'];

subj.fileName = ['subjects/',subj.ID,'_Position','_',date,'_',timeStr];
save(subj.fileName,'subj');