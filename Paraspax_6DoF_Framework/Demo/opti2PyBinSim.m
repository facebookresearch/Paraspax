%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function opti2PyBinSim(optiTimer);
%
% Function to stream optitrack data and send to pybinsim. For more details 
% on the Paraspax real-time frameweork, see [1].
%
% Output:
% -                   
%
% Input:        
% optiTimer         - timer object for Optitrack streaming
%
% INFO: The timer running the experiment can also change the osc command
% send in the nested function 'send2pybinsim' by changing the global variable
% binsim.currentSource. However, this seemed to be more safe than using a
% function to send osc commands which every timer can access, because then
% runtim errors could occure when both timers access the function the same
% time with different commands. The way it is implemented now, no locks are
% required, even though the software structure is a little harder to understand.
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

function opti2PyBinSim(optiTimer)
    
%Get global variables
global binsim

%Add required info to optiTimer
optiTimer.ErrorFcn = @(myTimerObj,thisEvent)warning('Something is not working...');
optiTimer.TimerFcn = @optiTimerFcn; %Run
optiTimer.StartFcn = @optiTimerStartFcn; %Init
optiTimer.StopFcn = @optiTimerStopFcn; %Cleanup possible

%Get grid data
cartGrid = binsim.cartGrid;
sphGrid = binsim.sphGrid;
%Convert spherical grid
sphGrid(:,2) = 90-sphGrid(:,2);
sphGridRad = sphGrid * pi / 180;
[xSphGrid,ySphGrid,zSphGrid] = sph2cart(sphGridRad(:,1),sphGridRad(:,2),1);

%Init variables required later
natnetclient = connect_to_optitrack();
%get information on the tracked objects:
model = natnetclient.getModelDescription;
%set listener information (and Optitrack tracking model names):
listeners = 1:model.RigidBodyCount;
%build cell array of available names:
names = cell(1,numel(listeners));
for listener = 1:numel(listeners)
    names{listener} = ['RigidBody ',sprintf('%02d',listeners(listener))];
end
%Pre-Define listener_data struct
for listener = 1:numel(listeners)
    listener_data(listener).name = [];

    %match to the list of names (find the index):
    listener_data(listener).name_idx = [];

    %get quaternion for listeners:
    listener_data(listener).quaternion = [];

    %get listener positions:
    listener_data(listener).position = [];
end
%Write everything in init struct and add to UserData of optiTimer
init.natnetclient = natnetclient;
init.model = model;
init.listeners = listeners;
init.names = names;
init.listener_data = listener_data;
%Init further variables
az = nan;
el = nan;
rotated_pos = nan(1,3);
relativeAngles = nan(2,1);
relativeAnglesRad = nan(2,1);
lpPar = nan(1,3);
xRelative = nan;
yRelative = nan;
zRelative = nan;
euDistance = nan;
angleID = nan;
posID = nan;
%Define id for headphone (here always one)
idListener = 1;

%Start timer
start(optiTimer);

%% Timer functions

%Initialize connection to Optitrack
function optiTimerStartFcn(~,~)

    disp('Streaming Optitrack Data...');

end

function optiTimerFcn(~,~)
    
    %---- GET TRACKING DATA ----

    %Grab new data:
    data = init.natnetclient.getFrame;
    
    idSource = binsim.currentSource+2; %Global variable which gets updated by midi device or experimental control, +2 because binsim starts with 0 and listene 1 here = headphones
    for listener = [idListener,idSource] %Optimization - Get just data of current listener-source condition. 1 is always headphone

        %get rigid body name:
        init.listener_data(listener).name = init.model.RigidBody(listener).Name;
                
        %match to the list of names (find the index):
        init.listener_data(listener).name_idx = find(contains(init.names,init.listener_data(listener).name));

        %get quaternion for listeners:
        init.listener_data(listener).quaternion = [data.RigidBody(listener).qx,...
            data.RigidBody(listener).qy,data.RigidBody(listener).qz,data.RigidBody(listener).qw];
        
        %get listener positions:
        init.listener_data(listener).position = [data.RigidBody(listener).x,data.RigidBody(listener).y,data.RigidBody(listener).z];

    end
    
    %establish relative position of other listeners:
    relative_pos = init.listener_data(idSource).position - init.listener_data(idListener).position;

    %re-order quaternion:
    listen_quat = quaternion([init.listener_data(idListener).quaternion(4),...
        init.listener_data(idListener).quaternion(1:3)]);
    
    %rotate the cartesian coordinates
    rotated_pos = listen_quat.rotateframe(relative_pos);

    %convert to az/el:
    [az,el] = cart2sph(rotated_pos(:,3),rotated_pos(:,1),rotated_pos(:,2));

    %convert azimuth to -degrees and correct for Matlab coordinate sys:
    relativeAngles(1) = mod((az'*180/pi),360)-180;

    %convert elevation to degrees and rotate up (positive angles only):
    relativeAngles(2) = mod((el'*180/pi)+180,360)-180;

    %---- CHANGE TO PARASPAX COORDINATE SYSTEM ----
    
    %Optitrack is (X,Y,Z), as definied in "Angle" room with -Z along
    %the room length (towards the loudspeakers), +X along the room width
    %(towards the windows), Y going along height. This convention sadly
    %depends on how the OptiTrack system was calibrated in the room, so
    %this has to be adjuste correctly all the time! Some typical use
    %cases could be implemented though!

    %Listener position for paraspax
    lpPar(1) = -init.listener_data(idListener).position(3);
    lpPar(2) = -init.listener_data(idListener).position(1);
    lpPar(3) = init.listener_data(idListener).position(2);

    %---- RELATIVE ANGLES ----
    %Adjust azimuth to 0-360
    relativeAngles(relativeAngles(1,:)<0) = relativeAngles(relativeAngles(1,:)<0)+360;
    %Convert to radiant
    relativeAnglesRad = relativeAngles * pi / 180;
    %Convert to cartesian coordinates
    [xRelative,yRelative,zRelative] = sph2cart(relativeAnglesRad(1,:),relativeAnglesRad(2,:),1);
    %Calculate Euclidean distance
    euDistance = sqrt(sum(([xSphGrid,ySphGrid,zSphGrid] - [xRelative,yRelative,zRelative]) .^ 2,2));
    %Get ID with smallest distance
    [~,angleID] = min(euDistance);

    %---- POSITION ----
    %Find nearest position for grid with Euclidean distance (without height)
    %euDistance = sqrt(sum(([lpPar(1),lpPar(2),lpPar(3)] - [cartGrid(:,1),cartGrid(:,2),cartGrid(:,3)]) .^ 2,2))
    euDistance = sqrt(sum(([lpPar(1),lpPar(2)] - [cartGrid(:,1),cartGrid(:,2)]) .^ 2,2));
    %Get ID with smallest distance (starts at 0)
    [~,posID] = min(euDistance); posID = posID-1;

    %---- OSC MESSAGE ----
    %Channel, AZ, EL, POSITION, SOURCE
    %CHANGES EARLY FILTERS - DIFFERENT COMMAND FOR LATE REVERB FILTERS (not changed in this application)
    %90-sphGrid(angleID,2) because Paraspax works with colatitude
    oscsend(binsim.u,'/pyBinSimFilter','iiiiiiiiii',binsim.channelID,sphGrid(angleID,1),90-sphGrid(angleID,2),posID,binsim.currentSource,0,0,0,0,0);
    
end

function optiTimerStopFcn(~,~)

    disp('Opti Timer Stopped...');

end

%Function to set the late reverb filter. 
%Not gonne be changed in this application, but with the ID, it could be
%adjusted in general independent to the source. In this application, it is
%coupled to binsim.currentSource. ReverbFilters must all have the same length!
function setReverbFilter()
  
    oscsend(binsim.u,'/pyBinSimLateReverbFilter','iiiiiiiiii',binsim.channelID,0,0,binsim.reverbFilterID,0,0,0,0,0,0);
   
end

end