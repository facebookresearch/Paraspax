%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function [equiangleGrid, nPoints] = getEquiangleGrid(azResolution,elResolution,excludePoles,truncAz,truncEl)
%
% This function can be used to calculate a equiangular grid of any type for
% spherical BRIR synthesis. The function creates a grid based on the 
% resolution in azimuth and elevation with equidistance between the grid
% points. Moreover, the grid always starts at 90 degree elevation
% (SH coordinates / colatitude), thus supporting frontal head orientation by default.
%
% Output:
% gaussGrid         - Array with corresponding gaussGrid in SH
%                     coordinates
%                     [azDEG, elDEG]  
% nPoints           - Total umber of grid points
%
% Input: 
% azResolution      - Azimuth resolution
%                     azResolution = 1 --> Steps of 1° in the horizontal
%                     plane from 0 - 359deg --> 360deg values per elevation
%                     step
%                     azResolution = 2 --> Steps of 2deg from 0 - 358deg...
% elResolution      - Elevation resolution
%                     elResolution = 10 --> Steps of 10deg in the median plane
%                     starting from 90deg up to 0deg and down to 180deg --> 19
%                     elevation values 
%                     elResolution = 30 --> Steps of 30deg --> 7 values
% excludePoles      - True/false - exclude values at poles (0deg and 180deg) 
%                     --> number of elevation values recuded by two if set
%                     true
%                     Default: true  
% truncAz           - +/- truncation value in the horizontal plane. For example 
%                     truncAz = 30 --> azRange from 30deg to 330deg in SH coordinates
%                     (-30deg to +30deg)
%                     Default: 0 - no truncation
%                     Use [a b] vector for non-symmetrical truncation, for
%                     example [60 330] would truncate values to this range
% truncEl           - +/- truncation value in the median plane. For example 
%                     truncEl = 30 --> elRange from 60deg to 120deg in SH coordinates
%                     (-30deg to +30deg)
%                     Use [a b] vector for non-symmetrical truncation, for
%                     example [30 120] would truncate values to this range
%                     Default: 0 - no truncation
%
% Dependencies: -
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.
%
% Copyright (c) Facebook, Inc. and its affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

function [equiangleGrid, nPoints] = getEquiangleGrid(azResolution,elResolution,excludePoles,truncAz,truncEl)

if nargin < 3 || isempty(excludePoles)
    excludePoles = true;
end

if nargin < 4 || isempty(truncAz)
    truncAz = 0;
end

if nargin < 5 || isempty(truncEl)
    truncEl = 0;
end

if azResolution < 1 || elResolution < 1
    error('The smallest possible resolution is 1!');
end

if (floor(azResolution) ~= azResolution) || (floor(elResolution) ~= elResolution)
    error('Resolution can only be defined in integers!');
end

if truncAz(1) < 0 || truncEl(1) < 0
    error('Sorry, but truncation values must be > 0!');
end

if size(truncAz,2) == 2
    if truncAz(2) < 0
        error('Sorry, but truncation values must be > 0!');
    end
end

if size(truncEl,2) == 2
    if truncEl(2) < 0
        error('Sorry, but truncation values must be > 0!');
    end
end

%%
%Get values according to resolution
%Azimuth values
azValues = 0:azResolution:359;
azValues = azValues';

if truncAz(1) > 0
    if size(truncAz,2) == 1 %Symmetrical truncation values
        azValues = [azValues(azValues<=truncAz);azValues(azValues>=360-truncAz)];
    end
    
    if size(truncAz,2) == 2 %Non-symmetrical truncation values
        azValues = [azValues(azValues<=truncAz(1));azValues(azValues>=truncAz(2))];
    end
end

%Elevation values
elValues1 = 90:elResolution:180;
elValues2 = 90:-elResolution:0;
elValues = [flip(elValues2),elValues1(2:end)];
elValues = elValues';

if excludePoles
    elValues = elValues(2:end-1);
end

if truncEl(1) > 0
    if size(truncEl,2) == 1 %Symmetrical truncation values
        elValues = (elValues(elValues<=90+truncEl));
        elValues = (elValues(elValues>=90-truncEl));
    end
    
    if size(truncEl,2) == 2 %Non-symmetrical truncation values
        elValues = elValues(elValues>=truncEl(1));
        elValues = elValues(elValues<=truncEl(2));
    end
end
   
%Combine to final grid in DEG
grid1 = repmat(azValues,[1,size(elValues)]);
grid1 = grid1';
grid1 = grid1(:);
grid2 = repmat(elValues,[size(azValues),1]);
equiangleGrid = [grid1,grid2];

%Get nPoints
nPoints = size(equiangleGrid,1);

end