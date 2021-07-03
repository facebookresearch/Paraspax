% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function pyBinSimMidiControl(midiTimer)
%
% Function to run and control experiment with Numark Orbit MIDI device. For 
% more details on the Paraspax real-time framework, see [1].
%
% Output:
% -     
%
% Input:        
% midiTimer         - timer object for MIDI streaming
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


function pyBinSimMidiControl(midiTimer)
    
%Get global struct settingsHandle
global sh;
global binsim;
 
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

%Chose nice test signal as first source
sh.testSignalID = find(strcmp(sh.testSigWavNames,'Singing_Male_02.wav')==1);

%Stop pybinsim, just to be sure
oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','False');

%Start timer
start(midiTimer);

%% Timer functions
function midiTimerStartFcn(~,~)

    midireceive(device);midireceive(device);midireceive(device);%Clear buffer
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

    if (sum(midiChannelID == 10)/size(msg,1)) == 1
        disp('Source 1');
        %Chose sourceID in global handle. This is gonne change the OSC
        %message send to pybinsim
        binsim.currentSource = binsim.sourceIDs(1);
        %Set reverbFilter accordingly
        oscsend(binsim.u,'/pyBinSimLateReverbFilter','iiiiiiiiii',binsim.channelID,0,0,binsim.currentSource,0,0,0,0,0,0);
        %Play again
        playAgain(); 
    end

    if (sum(midiChannelID == 11)/size(msg,1)) == 1
        disp('Source 2');
        binsim.currentSource = binsim.sourceIDs(2);
        %Set reverbFilter accordingly
        oscsend(binsim.u,'/pyBinSimLateReverbFilter','iiiiiiiiii',binsim.channelID,0,0,binsim.currentSource,0,0,0,0,0,0);
        %Play again
        playAgain(); 
    end

    if (sum(midiChannelID == 12)/size(msg,1)) == 1
        disp('Source 3');
        binsim.currentSource = binsim.sourceIDs(3);
        %Set reverbFilter accordingly
        oscsend(binsim.u,'/pyBinSimLateReverbFilter','iiiiiiiiii',binsim.channelID,0,0,binsim.currentSource,0,0,0,0,0,0);
        %Play again
        playAgain(); 
    end

    if (sum(midiChannelID == 13)/size(msg,1)) == 1
        disp('Source 4');
        binsim.currentSource = binsim.sourceIDs(4);
        %Set reverbFilter accordingly
        oscsend(binsim.u,'/pyBinSimLateReverbFilter','iiiiiiiiii',binsim.channelID,0,0,binsim.currentSource,0,0,0,0,0,0);
        %Play again
        playAgain(); 
    end

    if (sum(midiChannelID == 14)/size(msg,1)) == 1
        disp('Real / Virtual Toggle');
        
        toggleRealVirtual();
    end

    if (sum(midiChannelID == 15)/size(msg,1)) == 1
        disp('Random Test Signal');
       
        randomTestSignal();
        
    end
    
    if (sum(midiChannelID == 16)/size(msg,1)) == 1
        disp('Play again');
       
        playAgain();
        
    end

    if (sum(midiChannelID == 9)/size(msg,1)) == 1
        disp('Stop Playback');
       
        stopPlayback();
        
    end

end
end

function midiTimerStopFcn(~,~)
    disp('MIDI Timer Stopped');
end

function toggleRealVirtual()

    if sh.rv == 0 %If currently real source
        
        %Get analog output channel
        outputChannel = binsim.currentSource + 1;
        
        %Stop if current object played
        stop(audioPlayerObj);
        audioPlayerObj = audioplayer(zeros(16,2),sh.io.fs,sh.io.nBits,sh.io.outputIdMatrix(outputChannel));
        play(audioPlayerObj);  
        
        %Pause pybinsim
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True')
        %Fill with zeros
        oscsend(binsim.u,'/pyBinSimFile','s',zeroBufferCommand);
        
        %Set to virtual source
        sh.rv = 1;
               
        %Load current audiofile into pybinsim
        testSigWavName = sh.testSigWavNames{sh.testSignalID};
        oscCommand = ['signal/',testSigWavName];
        
        %Play stimulus
        oscsend(binsim.u,'/pyBinSimFile','s',oscCommand);
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','False')
        
    else %If currently virtual source
       
        %Pause pybinsim
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True')
        %Fill with zeros
        oscsend(binsim.u,'/pyBinSimFile','s',zeroBufferCommand);
        
        %Set to real source
        sh.rv = 0;
        %Get analog output channel
        outputChannel = binsim.currentSource + 1;
        %Get test signal
        testSignal = sh.testSig.(sh.testSigFieldNames{sh.testSignalID});
        
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
          
    end
    
end

function randomTestSignal()
    
    if sh.rv == 0 %Real source
        
        %Stop if current object played
        stop(audioPlayerObj);
        
        %Get random signal from pool
        sh.testSignalID = randi([1,length(sh.testSigWavNames)]);
        %Get analog output channel
        outputChannel = binsim.currentSource + 1;
        %Get test signal
        testSignal = sh.testSig.(sh.testSigFieldNames{sh.testSignalID});
        
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

    else %Virtual source
        
        %Pause pybinsim
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True')
        %Fill with zeros
        oscsend(binsim.u,'/pyBinSimFile','s',zeroBufferCommand);
        
        %Get random signal from pool
        sh.testSignalID = randi([1,length(sh.testSigWavNames)]);
        
        %Load audiofile into pybinsim
        testSigWavName = sh.testSigWavNames{sh.testSignalID};
        oscCommand = ['signal/',testSigWavName];
        
        %Play stimulus
        oscsend(binsim.u,'/pyBinSimFile','s',oscCommand);
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','False')
    
    end
end

function playAgain()
    
    if sh.rv == 0 %If currently real source
        
        %Stop if current object played
        stop(audioPlayerObj);
        
        %Get analog output channel
        outputChannel = binsim.currentSource + 1;
        %Get test signal
        testSignal = sh.testSig.(sh.testSigFieldNames{sh.testSignalID});
        
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
        
    else %If currently virtual source

        %Pause pybinsim
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True')
        %Fill with zeros
        oscsend(binsim.u,'/pyBinSimFile','s',zeroBufferCommand);
        
        %Load current audiofile into pybinsim
        testSigWavName = sh.testSigWavNames{sh.testSignalID};
        oscCommand = ['signal/',testSigWavName];
        
        %Play stimulus
        oscsend(binsim.u,'/pyBinSimFile','s',oscCommand);
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','False')
        
    end
    
end

function stopPlayback()
    
    if sh.rv == 0 %If currently real source
        
        %Stop if current object played
        stop(audioPlayerObj);
                
    else %If currently virtual source
        
        %Pause pybinsim
        oscsend(binsim.u,'/pyBinSimPauseAudioPlayback','s','True')
        %Fill with zeros
        oscsend(binsim.u,'/pyBinSimFile','s',zeroBufferCommand);
        
    end
     
end
end