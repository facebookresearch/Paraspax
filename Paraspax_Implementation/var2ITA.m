% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function ITAobj = var2ITA(var,fs)
%
% Function to convert a variable (usually a RIR) in workspace 
% to an ITA object for further room acoustical analysis with the ITA
% toolbox
%
% Output:
% ITAobj                - ITA object
%
% Input:        
% varname               - Variable in workspace
% fs                    - Sample rate of variable / RIR
%
% Dependencies: ITA toolbox
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function ITAobj = var2ITA(var,fs)

ITAobj = itaAudio;
ITAobj.samplingRate = fs;
ITAobj.timeData = var;

end