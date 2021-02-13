%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function psx = extrapolate(psx)
%
% Function to extrapolate parameters of spatialized RIR to new listener
% position. Extrapolation defined by X Y Z vector leading to a relative 
% shift of receiver in relation to source.
%
% Output:
% psx                   - psx struct with extrapolated parameters in 
%                         psx.extrap field
%
% Input:        
% psx                   - psx struct with required fields
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

function psx = extrapolate(psx)

%Get various variables for further processing
fs  = psx.fs;
c   = psx.c; 
lps = psx.extrap.lps;

%Head orientation of listener, defined as az/el.
headOrientation = [0,0];

%% Perform extrapolation - listener translation

%Get distance of direct sound and reflections in meter
rDirect = psx.geo.sourceDistance;
rRef = psx.spat.refListSpat.toa / fs * c; 

%Get az/el of direct sound and reflections and transform to radiant
%Dont use SH angles!
azElDirect = psx.geo.sourceDirection;
azElDirect = azElDirect * pi / 180;
azElRef = psx.spat.refListSpat.selectAzEl;
azElRef = azElRef * pi / 180;

%Combine everything in one array
azElR = [[azElDirect;azElRef],[rDirect;rRef]];

%Transform to cartesian coordinates
%This are relativ coordinates between source and receiver, not absolute
%coordinates of sourcer/receiver!
[x,y,z] = sph2cart(azElR(:,1),azElR(:,2),azElR(:,3));

%Just a small check if transformation was correct
if psx.geo.roomDimensions
    if sum(round([psx.geo.srcPos-psx.geo.recPos]-[x(1) y(1) z(1)],3)) ~= 0
        warning('Transformation of spherical to cartesian coordinates might be incorrect');
    end
end

%Extrapolate listener position
%Negative sign because shift is applied to listener and X Y Z are related
%to source
xExtrap = x-lps(1);
yExtrap = y-lps(2);
zExtrap = z-lps(3);

%Transform back to spherical coordinates
[azElRextrap(:,1),azElRextrap(:,2),azElRextrap(:,3)] = cart2sph(xExtrap,yExtrap,zExtrap);
azElRextrap(:,1) = azElRextrap(:,1) / pi * 180;
azElRextrap(:,1) = azElRextrap(:,1) + 360;
azElRextrap(:,1) = mod(azElRextrap(:,1),360);
azElRextrap(:,2) = azElRextrap(:,2) / pi * 180;

%Calculate distance difference
%Negative sign of rShift means reduced distance.
rExtrap = azElRextrap(:,3) - azElR(:,3);

%Calculate amplitude factors
ampFactorsExtrap = azElR(:,3) ./ azElRextrap(:,3);

%Get extrapolated TOA vector in samples
rExtrapSamples = round(azElRextrap(:,3) * fs / c);

%Get new start and end sample points 
SE_Direct_Extrap = [rExtrapSamples(1)-round((psx.refDet.directWin(1)*fs)), rExtrapSamples(1)+round((psx.refDet.directWin(2)*fs))];
SE_Ref_Extrap = [rExtrapSamples(2:end)-round((psx.refDet.refWin(1)*fs)), rExtrapSamples(2:end)+round((psx.refDet.refWin(2)*fs))];

%Make some checks. Could be handled in various ways...
if SE_Direct_Extrap(1) < 1
    warning('Direct sound window has negative samples. Add safety margin at beginning of RIR for BRIR synthesis!')
end
if SE_Ref_Extrap(:,1) < 1
    warning('Reflection window has negative samples. Add safety margin at beginning of RIR for BRIR synthesis!')
end
if rExtrapSamples(2:end) < rExtrapSamples(1)
    warning('TOA of reflection(s) shorter than TOA of direct sound!')
end

%Apply rotation to sources according to head orientation
if sum(abs(headOrientation)) ~= 0
    [azElRextrap(:,1),azElRextrap(:,2)] = AKroomSimulationRotation(azElRextrap(:,1),azElRextrap(:,2),headOrientation(1),headOrientation(2));
end
    
%% Write results in struct

if psx.geo.roomDimensions
    psx.extrap.recPosExtrap = psx.geo.recPos + lps;
end
psx.extrap.refListSpat.tDirect = rExtrapSamples(1);
psx.extrap.refListSpat.rExtrapDirect = rExtrap(1);
psx.extrap.refListSpat.ampFactorDirect = ampFactorsExtrap(1);
psx.extrap.refListSpat.AzEl_Direct = azElRextrap(1,1:2);
psx.extrap.refListSpat.AzEl_Direct_SH = psx.extrap.refListSpat.AzEl_Direct;
psx.extrap.refListSpat.AzEl_Direct_SH(:,2) = 90-psx.extrap.refListSpat.AzEl_Direct_SH(:,2);
psx.extrap.refListSpat.SE_Direct = SE_Direct_Extrap;
psx.extrap.refListSpat.SrcAzEl_Direct = [azElRextrap(1,1), -azElRextrap(1,2)];%From the point of view of the source
psx.extrap.refListSpat.SrcAzEl_Direct_SH = psx.extrap.refListSpat.SrcAzEl_Direct;
psx.extrap.refListSpat.SrcAzEl_Direct_SH(:,2) = 90 - psx.extrap.refListSpat.SrcAzEl_Direct_SH(:,2);
psx.extrap.refListSpat.srcView = [mod(180+azElRextrap(1,1), 360) -azElRextrap(1,2)];
psx.extrap.refListSpat.toa = rExtrapSamples(2:end);
psx.extrap.refListSpat.rExtrap = rExtrap(2:end);
psx.extrap.refListSpat.ampFactor = ampFactorsExtrap(2:end);
psx.extrap.refListSpat.rExtrap = rExtrap(2:end);
psx.extrap.refListSpat.selectAzEl = azElRextrap(2:end,1:2);
psx.extrap.refListSpat.selectAzEl_SH = psx.extrap.refListSpat.selectAzEl;
psx.extrap.refListSpat.selectAzEl_SH(:,2) = 90-psx.extrap.refListSpat.selectAzEl(:,2); 
psx.extrap.refListSpat.selectSE = SE_Ref_Extrap;

end