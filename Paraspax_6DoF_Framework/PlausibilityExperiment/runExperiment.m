% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function function runExperiment(midiTimer)
%
% Function to run and control the 6-DoF plausibility experiment with MIDI device.
% Gets called from Paraspax_PlausibilityExperiemnt script. For more details on the 
% Paraspax real-time framework and the plausibility study, see [1].
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

function runExperiment(midiTimer)
    
%Get global struct settingsHandle
global sh;
global binsim;
global subj;
 
%Define rest of midiTimer
midiTimer.ErrorFcn = @(myTimerObj,thisEvent)warning('Error');
midiTimer.TimerFcn = @midiTimerFcn; %Run
midiTimer.StartFcn = @midiTimerStartFcn; %Init
midiTimer.StopFcn = @midiTimerStopFcn; %Cleanup
    
%Get global variables in nested function
device = midiTimer.UserData.device;
midiChannelID = [];
msg = midimsg(0);

%Init audioplayerobject
audioPlayerObj = audioplayer(zeros(64,2),sh.io.fs,sh.io.nBits,sh.io.outputIdMatrix(1));

%Zero buffer for pybinsim
zeroBuffer = zeros(16,1); audiowrite('zeroBuffer.wav',zeroBuffer,sh.io.fs);
zeroBufferPath = pwd;
zeroBufferCommand = [zeroBufferPath,'\zeroBuffer.wav'];

%Stop pybinsim, just to be sure
oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True');

%Boolean isPlaying
virtualIsPlaying = 0;

%Start timer
start(midiTimer);

%% Timer functions
function midiTimerStartFcn(~,~)
    
    midireceive(device);midireceive(device);midireceive(device);%Clear buffer
    msg = midimsg(0);
    disp('Streaming MIDI Data...'); 
    
end

function midiTimerFcn(~,~)
    
    msg = midimsg(0);%Reset msg
    midiChannelID = []; %Reset midiChannelID
    msg = midireceive(device);
    
    if size(msg,1) > 0 %Not empty
       
        %Channels
        for kk = 1:size(msg,1)
            midiChannelID(kk) = msg(kk).Channel;
        end
               
        %Check if only channel 1/2 chosen to be sure
        if (sum(midiChannelID == 1)/size(msg,1)) == 1
            disp('Left Pads - Blue');  
            
            if virtualIsPlaying
                oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True');
                virtualIsPlaying = 0;
            end    
            if isplaying(audioPlayerObj)
               stop(audioPlayerObj); 
            end
            
            %Test starts with any keypress. Dont save data
            if subj.trial == 0 
                disp('Test started');
                runTrial(); 
                
            else %If not test start, get response and save data
                          
                resp = 0; %0 - Left / Real
                subj.resp(subj.trial) = resp;
                
                %Check if answer correct
                if resp == subj.rvsequ(subj.trial) %If correct
                    subj.respCorrect(subj.trial) = 1;
                else
                    subj.respCorrect(subj.trial) = 0;
                end
                
                %Save before next trial
                save(subj.fileName,'subj');
                
                %Run next trial if trial left
                if subj.trial < subj.nTrials
                    runTrial();
                else %If no trials left, end experiment
                    endExperiment();
                end
            end  
        end
        
        if (sum(midiChannelID == 2)/size(msg,1)) == 1
            disp('Right Pads - Green');
            
            if virtualIsPlaying
                oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True');
                virtualIsPlaying = 0;
            end    
            if isplaying(audioPlayerObj)
               stop(audioPlayerObj); 
            end
            
            %Test starts with any keypress. Dont save data
            if subj.trial == 0 
                disp('Test started');
                runTrial();
                
            else %If not test start, get response and save data
                
                resp = 1; %1 - Right / Virtual
                subj.resp(subj.trial) = resp; 
                
                %Check if answer correct
                if resp == subj.rvsequ(subj.trial) %If correct
                    subj.respCorrect(subj.trial) = 1;
                else
                    subj.respCorrect(subj.trial) = 0;
                end
                
                %Save before next trial
                save(subj.fileName,'subj');
                
                %Run next trial if trial left
                if subj.trial < subj.nTrials
                    runTrial();
                else %If no trials left, end experiment
                    endExperiment();
                end
            end
        end
        
    end
end

function midiTimerStopFcn(~,~)
    disp('MIDI Timer Stopped');
end

function runTrial()
    
    %Stop midiTimer so no input can be done while trial is played
    %stop(midiTimer); 
    msg = midimsg(0);%Reset msg
    
    %Stop pybinsim just to be sure
    oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True');
    %Fill with zeros
    oscsend(binsim.u,'/pyBinSimFile','s',zeroBufferCommand);
    
    %Clear audioplayer object
    stop(audioPlayerObj);
    for kk = 1:sh.nSources
        audioPlayerObj = audioplayer(zeros(16,2),sh.io.fs,sh.io.nBits,sh.io.outputIdMatrix(kk));
        play(audioPlayerObj);
    end
    
    %Increase trial number
    subj.trial = subj.trial + 1;
    
    %Show trial number in console
    fprintf('Trial %d\n',subj.trial);
    
    %Get condition
    %(1) - Check if real or virtual source
    if subj.rvsequ(subj.trial) == 0 %Real source condition
        
        %Clear audioplayer object
        %audioPlayerObj = [];
        
        %Get analog output channel for condition
        outputChannel = subj.sourceSequ(subj.trial) + 1;
        %Get test signal for condition
        testSignalID = subj.testSigSequ(subj.trial);
        testSignal = sh.testSig.(sh.testSigFieldNames{testSignalID});
        %tTestSignal = length(testSignal);
                
        %Audioplayer always plays both channels. Add zeros...
        if outputChannel == 1 || outputChannel == 3
            testSignal = [testSignal,zeros(length(testSignal),1)];
        end
        if outputChannel == 2 || outputChannel == 4
            testSignal = [zeros(length(testSignal),1),testSignal];
        end
        
        %Create audioPlayerObject and play
        audioPlayerObj = audioplayer(testSignal,sh.io.fs,sh.io.nBits,sh.io.outputIdMatrix(outputChannel));
        play(audioPlayerObj)
        
        %Start midiTimer again after playback
        %midiTimer.StartDelay = round(tTestSignal/sh.io.fs+0.5,2);
        %start(midiTimer);
                
    end
    
    if subj.rvsequ(subj.trial) == 1 %Virtual source condition
        
        %Get ID for virtual current source and update global variable send
        %to pybinsim in opti2PyBinSim
        binsim.currentSource = subj.sourceSequ(subj.trial); %Starts at 0
        %Set reverbFilter accordingly
        oscsend(binsim.u,'/pyBinSimLateReverbFilter','iiiiiiiiii',binsim.channelID,0,0,binsim.currentSource,0,0,0,0,0,0);
        
        %Load audiofile into pybinsim according to condition
        testSignalID = subj.testSigSequ(subj.trial);
        %tTestSignal = length(sh.testSig.(sh.testSigFieldNames{testSignalID}));
        testSigWavName = sh.testSigWavNames{testSignalID};
        oscCommand = ['signal/',testSigWavName];
        
        %Play stimulus
        virtualIsPlaying = 1;
        oscsend(binsim.u,'/pyBinSimFile','s',oscCommand);
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','False')
        
        %Start midiTimer again after playback
        %midiTimer.StartDelay = round(tTestSignal/sh.io.fs+0.5,2);
        %start(midiTimer);
        
    end
    
    %For evaluation
    trial = subj.trial
    rvsequ = subj.rvsequ(subj.trial)
    testSigSequ = subj.testSigSequ(subj.trial)
    sourceSequ = subj.sourceSequ(subj.trial)
        
    %Another pause just to be sure the midi timer does not start to early
    %pause(0.5);
      
    %Update global amount of trials
    %subj.trial = currentTrial;

end

function endExperiment()

    %Set to completed
    subj.expComp = 1;
    
    %Stop midiTimer
    stop(midiTimer);
    
    %Close log
    fclose(subj.logtxt);
    
    %Clear audioplayer object
    for kk = 1:sh.nSources
        audioPlayerObj = audioplayer(zeros(16,2),sh.io.fs,sh.io.nBits,sh.io.outputIdMatrix(kk));
        play(audioPlayerObj);
    end
    
    disp('Experiment complete...');
    %Play message indicating that experiment is complete
    expComplete = audioread('ExpComplete.wav');
    %AKio(expComplete,2,[],sh.io.fs,sh.io.playDev)
    
    outputChannel = 2;
    %Audioplayer always plays both channels. Add zeros...
    if outputChannel == 1 || outputChannel == 3
        expComplete = [expComplete,zeros(length(expComplete),1)];
    end
    if outputChannel == 2 || outputChannel == 4
        expComplete = [zeros(length(expComplete),1),expComplete];
    end

    %Create audioPlayerObject and play
    audioPlayerObj = audioplayer(expComplete,sh.io.fs,sh.io.nBits,sh.io.outputIdMatrix(outputChannel));
    play(audioPlayerObj)
          
end
end