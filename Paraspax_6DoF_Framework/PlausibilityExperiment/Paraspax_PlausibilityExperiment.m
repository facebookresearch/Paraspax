%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% script Paraspax_PlausibilityExperiment
%
% Script to run the 6-DoF plausibility experiment for the Paraspax project.
% For more details on the Paraspax real-time frameweork and the Plausibility 
% study, see [1].
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

%% START - (1) Open "Calibration Paraspax Exceptional (MeanErr 0.445 mm)2020-01-27 11" to start Optitrack
winopen('D:\Paraspax\Calibration Paraspax Exceptional (MeanErr 0.445 mm) 2020-01-27 11.cal')

%% START - (2) Run batch script "startPyBinSim" to start rendering engine and load filters (this takes about 30-60 minutes!)
winopen('D:\Paraspax\PyBinSimProjects\startPyBinSim_Angel_SDM10_Az4EL10_XYZ025.bat');

%% START - (3) ADD SUBJECT

global subj %Global variable so results can be saved there. Dont change

subj.ageRange = 2; %Age-range 1: < 25 || 2: 25 - 34 || 3: 35 - 44 || 4: 45 - 55 || 5: >55
subj.ID = 'ADRA_PSX_TEST'; 
subj.run = 0; 

%% If software crashed, load subj struct here

%global subj;
%subj = importdata('subj.mat'); subj.trial = subj.trial - 1; %Go one trial back just to be sure

%% Initialization of global variable - Dont change if you are not sure!

global sh; %Settings handle experiments
global binsim; %PyBinSim handle

%Test settings
sh.nTestSignals = 30;
sh.nSources = 4;
sh.pybinPath = 'D:\Paraspax\PyBinSimProjects\Angel_SDM10_Az4El10_XYZ025';
sh.expPath = pwd;

%PyBinSim settings
binsim.u = udp('127.0.0.1', 10000); fopen(binsim.u); %Open udp port directly
binsim.sourceIDs = [1:sh.nSources]-1;
binsim.currentSource = binsim.sourceIDs(2);
binsim.channelID = 0; %Always 0 for this experiment
binsim.reverbFilterID = 0; %Always 0 for this experiment
binsim.cartGrid = importdata([sh.pybinPath,'/grids/cartGrid.mat']);
binsim.sphGrid = importdata([sh.pybinPath,'/grids/sphGrid.mat']);

%Audio settings for AudioPlayer
sh.io.fs = 48000;
sh.io.nBits = 24;
sh.io.audioInfo = audiodevinfo;
for outputChannels = 1:size(sh.io.audioInfo.output,2)
   %tmp1 = strcmp(sh.io.audioInfo.output(outputChannels).Name,'Speakers (RME Fireface UCX) (Windows DirectSound)'); 
   tmp1 = strcmp(sh.io.audioInfo.output(outputChannels).Name,'Speakers (2- RME Fireface UCX) (Windows DirectSound)'); %New fireface interface
   %tmp2 = strcmp(sh.io.audioInfo.output(outputChannels).Name,'Analog (3+4) (RME Fireface UCX) (Windows DirectSound)');
   tmp2 = strcmp(sh.io.audioInfo.output(outputChannels).Name,'Analog (3+4) (2- RME Fireface UCX) (Windows DirectSound)'); %New fireface interface
   if tmp1
       sh.io.idCH12 = sh.io.audioInfo.output(outputChannels).ID;
   end
    if tmp2
       sh.io.idCH34 = sh.io.audioInfo.output(outputChannels).ID;
   end
end
sh.io.outputIdMatrix = [sh.io.idCH12,sh.io.idCH12,sh.io.idCH34,sh.io.idCH34];
clear tmp1 tmp2 outputChannels;

%Generarte new sequence for subject if test not already started
if subj.run == 0
    genSubjCond(); %Works with global variable subj
end

%% Load test signals from pybinsim folder

%Save in struct if signals are not exactly the same size
cd(sh.pybinPath);
cd('signal');
fileList = dir('*.wav');
testSig = [];
for kk = 1:numel(fileList)
    
    if kk < 10
        fieldName = ['sig0',num2str(kk)];
    else
        fieldName = ['sig',num2str(kk)];
    end
   
    testSigWavNames{kk} = fileList(kk).name;
    
    testSig(1).(fieldName) = audioread(fileList(kk).name);
end
sh.testSig = testSig;
sh.testSigFieldNames = fieldnames(testSig);
sh.testSigWavNames = testSigWavNames;

%Go back to experiment folder
cd(sh.expPath);
clear testSig testSigWavNames fieldName fileList kk;

%% Prepare timer

%MIDI settings - Timer also used to run entire experiment
midiTimer = timer;
midiTimer.Name = 'midiTimer';
midiTimer.BusyMode = 'drop';
midiTimer.ExecutionMode = 'fixedRate';
midiTimer.StartDelay = 0;
midiTimer.Period = 0.06; %approx. 15Hz is more than enough for MIDI
midiTimer.UserData.device = mididevice('Input','Numark ORBIT'); %Define device, can be checked with function mididevinfo. Save in UserData so timer can access it

%Optitrack settings
optiTimer = timer; 
optiTimer.Name = 'optiTimer';
optiTimer.BusyMode = 'drop';   
optiTimer.ExecutionMode = 'fixedRate';
optiTimer.Period = 0.006; %approx. 150Hz for tracking data

%% Real time streaming Optitrack

opti2PyBinSim_Log(optiTimer)

%% Start experiment

runExperiment(midiTimer);

%% Run in the end to clear timer

stop(optiTimer); 
stop(midiTimer); 

delete(optiTimer);
delete(midiTimer);
