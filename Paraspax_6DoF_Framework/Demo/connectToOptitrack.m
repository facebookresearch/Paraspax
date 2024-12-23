% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function connectToOptitrack();
%
% Function to connect to the Natnet OptiTrack SDK
%
% Output:
% -                   
%
% Input:        
% -
%
% Dependencies: Natnetbib from OptiTrack, OptiTrack
%
% References:
% [1] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on One Monaural Room Impulse Response,” 
% J. Audio Eng. Soc., vol. 69, no. 7/8, pp. 557–575, 2021. ﻿
% https://doi.org/10.17743/jaes.2021.0009
%
% Code written 2019/2020 by JMA, Johannes M. Arend.


function natnetclient = connectToOptitrack()

% Connect to NatNet:
natnetclient = natnet;

% Connect the client to the server
% Adjust for the respective network
fprintf( 'Connecting to the server\n' )
natnetclient.HostIP = '127.0.0.1';
natnetclient.ClientIP = '127.0.0.1';
natnetclient.ConnectionType = 'Multicast';
natnetclient.connect;
if ( natnetclient.IsConnected == 0 )
    fprintf( 'Client failed to connect\n' )
    fprintf( '\tMake sure the host is connected to the network\n' )
    fprintf( '\tand that the host and client IP addresses are correct\n\n' )
    return
end

% Get the model description for the names.
modelDescription = natnetclient.getModelDescription;
if ( modelDescription.RigidBodyCount < 1 )
    return
end

end