% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% script TEST_ISM_MOVE
%
% Script to run the 6-DoF plausibility experiment for the Paraspax project.
% For more details on the Paraspax real-time frameweork and the Plausibility 
% study, see [1]. Description see PlausibilityExperiment_Docs.
%
% Dependencies: oscsend, global variable 'sh' and 'binsim' providing
% information about udp and cartesian/spherical grid, Natnetbib from
% OptiTrack, Paraspax real-time framework, PyBinSim, OptiTrack
%
% References:
% [1] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on 
% One Monaural Room Impulse Response,�? Submitted for publication, 2020.
%
% Code written 2019/2020 by JMA, Johannes M. Arend.


%% START - (1) Open "Calibration Paraspax Exceptional (MeanErr 0.445 mm)2020-01-27 11" to start Optitrack
winopen('D:\Paraspax\Calibration Paraspax Exceptional (MeanErr 0.445 mm) 2020-01-27 11.cal')

%% START - (2) Run batch script "startPyBinSim" to start rendering engine and load filters (this takes about 15-30 minutes!)
winopen('D:\Paraspax\PyBinSimProjects\startPyBinSim_RoomA_SDM10_Az4EL10_XYZ025.bat');

%% START - (3) ADD SUBJECT

%BE SURE THAT YOU ARE IN THE FOLDER HAVING THE SAME NAME AS THIS SCRIPT

global subj %Global variable so results can be saved there.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subj.ageRange = 2; %Age-range 1: < 25 || 2: 25 - 34 || 3: 35 - 44 || 4: 45 - 55 || 5: >55
subj.ID = 'ADRA_PSX_ISM_MOVE_ID'; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subj.run = 0; 

%% If software crashed, load subj struct here

%global subj;
%subj = importdata('subjects/ADRA_PSX_ISM_MOVE_XXX.mat'); subj.trial = subj.trial - 1; %Go one trial back just to be sure

%% START - (4) Initialization of global variable - Dont change if you are not sure!

global sh; %Settings handle experiments
global binsim; %PyBinSim handle

%Test settings
sh.nTestSignals = 30;
sh.nSources = 4;
sh.pybinPath = 'D:\Paraspax\PyBinSimProjects\RoomA_SDM10_Az4El10_XYZ025';
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
   tmp1 = strcmp(sh.io.audioInfo.output(outputChannels).Name,'Speakers (2- RME Fireface UCX) (Windows DirectSound)');
   tmp2 = strcmp(sh.io.audioInfo.output(outputChannels).Name,'Analog (3+4) (2- RME Fireface UCX) (Windows DirectSound)');
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

%%%% LOAD TEST SIGNALS %%%%
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

%%%% PREPARE TIMER %%%%
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

%%%% START OPTITRACK STREAMING %%%%
opti2PyBinSim_Log(optiTimer)
disp('Experiment ready to start...');

%% (5) - Start experiment 

%Befor start check if:
%(1) - PyBinSim project completely loaded
%(2) - All cameras in Optitrack Move software are working properly
%(3) - Audio interface, loudspeakers, transmitter, and receiver are on
%(4) - Totalmix still runs with snapshot "PyBinSim"
%(5) - Receiver has enough battery
%(6) - MIDI device is on and on Pad Bank 1
%(7) - Participant is in listening area or at position 1

runExperiment(midiTimer);

%% (6) - Run in the end to clear timer

stop(optiTimer); 
stop(midiTimer); 
delete(optiTimer);
delete(midiTimer);
