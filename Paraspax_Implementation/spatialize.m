% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function psx = spatialize(psx)
%
% Function to spatialize the monaural RIR by asigning DOAs to the detected 
% early reflections, using pseudo-randomized DOAs, an image source model 
% (if geometrical information available) or predetermined DOAs, e.g., from
% SDM measurements (for more details, see [1]).
%
% Output:
% psx                   - psx struct with spatialized reflections / 
%                         new reflection list in psx.ISM field
%
% Input:        
% psx                   - psx struct with required fields
%
% Dependencies: AKtools
%
% References:
% [1] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on 
% One Monaural Room Impulse Response,�? Submitted for publication, 2020.
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function psx = spatialize(psx)

%Get various variables for further processing
fs  = psx.fs;
refListDet = psx.refDet.refListDet;
nRef = psx.refDet.nRef;
tDirect = psx.par.mon.rirOnset_20dB; %TOA of direct sound in samples
directWin = round(psx.refDet.directWin*fs); %Direct sound window in samples
refWin = round(psx.refDet.refWin*fs); %Direct sound window in samples
eDirect = psx.refDet.eDirect; %Energy / RMS of direct sound
if ~isfield(psx.geo,'doa') %Probably not available in SDM applied
    sourceDirection = psx.geo.sourceDirection; %Direction [az el] of direct sound
end
%% Perform ISM simulation if required geometrical information is available and if intended

%If no DOA information but room dimensions available
if ~isfield(psx.geo,'doa') && psx.geo.roomDimensions
    
    %Write room simulation (rs) struct for ISM simulation
    rs.L = [psx.geo.roomLength, psx.geo.roomWidth, psx.geo.roomHeight];
    rs.f = [.063 .08 .1 .125 .16 .2 .25 .315 .4	.5 .63	.8 1 1.25 1.6 2	2.5 3.15 4 5 6.3 8]' * 1e3;
    %Wall absorption coefficients alpha. Set to 1, as we just need basic 
    %image sources without frequency dependent characteristics
    %Values for 6 surfaces (walls).
    rs.alpha = ones(length(rs.f),6);
    rs.recPos  = psx.geo.recPos;
    rs.recView = [0 0]; %Viewing direction of the receiver
    rs.recRot = false; %Head rotation of receiver. Rotation is gonne be applied later
    rs.srcCoordinates = 'carthesian';
    rs.srcPos = psx.geo.srcPos;
    rs.srcView = [];
    rs.ISMtruncation = {'N' psx.geo.ISMorder}; %ISM order
    rs.c = psx.c; %Speed of sound
    
    %Calculate srcView facing the receiver
    rs.srcView = [mod(180+sourceDirection(:,1), 360) -sourceDirection(:,2)];
    
    %Run ISM simulation
    [ISM.d, ISM.A, ISM.Raz, ISM.Rel, ISM.Saz, ISM.Sel, ISM.X, ISM.Xrel, ISM.N, ISM.wallLog] = AKism(rs.L, rs.srcPos, rs.srcView, rs.recPos, rs.recView, rs.alpha, rs.ISMtruncation, rs.c);
    
    %Save in spatialization struct
    psx.spat.ISM = ISM;
    psx.spat.ISM.rs = rs;
    psx.spat.ISM.comment = 'Plot room with AKroomSimulationPlot(rs) and results with AKroomSimulationPlot(rs, ISM, 1, N (N in apostrophe) , 0:2)'; 
    
    
    %-------------------------------------%
    %Combine ISM with detected reflections
    %Add direct sound, could be from ISM too, but lets say we just measured
    %distance and azimuth/elevation in the room
    refListSpat.tDirect = tDirect;
    refListSpat.eDirect = eDirect;
    refListSpat.ampFactorDirect = 1;
    refListSpat.AzEl_Direct = sourceDirection;
    refListSpat.AzEl_Direct_SH = refListSpat.AzEl_Direct;
    refListSpat.AzEl_Direct_SH(:,2) = 90-refListSpat.AzEl_Direct_SH(:,2);
    refListSpat.SE_Direct = [tDirect-directWin(1),tDirect+directWin(2)];
    refListSpat.SrcAzEl_Direct = [sourceDirection(1), -sourceDirection(2)];%From the point of view of the source
    refListSpat.SrcAzEl_Direct_SH = refListSpat.SrcAzEl_Direct;
    refListSpat.SrcAzEl_Direct_SH(:,2) = 90 - refListSpat.SrcAzEl_Direct_SH(:,2);
    refListSpat.srcView = rs.srcView;
    
    %Get 1st and 2nd order reflection TOA
    idN1 = nan(length(ISM.d));
    idN2 = nan(length(ISM.d));
    for kk = 1:length(ISM.d)
        if ISM.N(kk) == 1
            idN1(kk) = kk;
        end
        if ISM.N(kk) == 2
            idN2(kk) = kk;
        end
    end
    idN1 = idN1(~isnan(idN1));
    idN2 = idN2(~isnan(idN2));

    %d, dN1 and dN2 in samples
    d = round(ISM.d*fs/rs.c);
    dN1 = ISM.d(idN1); dN1 = round(dN1*fs/rs.c);
    dN2 = ISM.d(idN2); dN2 = round(dN2*fs/rs.c);

    %Save in struct
    refListSpat.toa = refListDet(1:nRef,1); %TOA/E from detection list
    refListSpat.e   = refListDet(1:nRef,2);
    refListSpat.ampFactor   = ones(nRef,1);
    
    %Check first order reflections first
    refListSpat.N1 = nan(length(refListSpat.toa),1);
    refListSpat.AzEl_N1 = nan(length(refListSpat.toa),2);
    refListSpat.SE_N1 = nan(length(refListSpat.toa),2);
    for kk = 1:length(dN1)
        td = refListSpat.toa-dN1(kk);
        [~,id] = min(abs(td));
        id2 = idN1(kk);
        refListSpat.N1(id) = id2;
        
        refListSpat.AzEl_N1(id,:)  = [ISM.Raz(id2),ISM.Rel(id2)];
        refListSpat.SE_N1(id,:)    = [refListDet(id,1)-refWin(1),refListDet(id,1)+refWin(2)];
    end
    
    %Then go through second order reflections
    refListSpat.N2 = nan(length(refListSpat.toa),1);
    refListSpat.AzEl_N2 = nan(length(refListSpat.toa),2);
    refListSpat.SE_N2 = nan(length(refListSpat.toa),2);
    for kk = 1:length(dN2)
        td = refListSpat.toa-dN2(kk);
        [~,id] = min(abs(td));
        id2 = idN2(kk);
        refListSpat.N2(id) = id2;

        refListSpat.AzEl_N2(id,:) = [ISM.Raz(id2),ISM.Rel(id2)];
        refListSpat.SE_N2(id,:)   = [refListDet(id,1)-refWin(1),refListDet(id,1)+refWin(2)];
    end
    
    %Last get the global minimum distance no matter which order, but still
    %sort according to order
    refListSpat.GN1 = nan(length(refListSpat.toa),1);
    refListSpat.AzEl_GN1 = nan(length(refListSpat.toa),2);
    refListSpat.SE_GN1 = nan(length(refListSpat.toa),2);
    refListSpat.GN2 = nan(length(refListSpat.toa),1);
    refListSpat.AzEl_GN2 = nan(length(refListSpat.toa),2);
    refListSpat.SE_GN2 = nan(length(refListSpat.toa),2);
    usedID = nan(length(refListSpat.toa),1);
    for kk = 1:length(refListSpat.toa)
        td = refListSpat.toa(kk)-d;
        [~,id] = min(abs(td));
        
        %Check if id was alread used, and assign next closest ID if it was
        %already used...
        zz = 1;
        while ~isempty(intersect(id,usedID))
            zz = zz+1;
            nc = unique(abs(td));
            ncID = nc(zz);
            id = find(abs(td) == ncID);
        end
        %Finally save ID not used before
        usedID(kk) = id;    
        
        %Save according to N
        if ISM.N(id) == 1
            refListSpat.GN1(kk) = id;
            refListSpat.AzEl_GN1(kk,:) = [ISM.Raz(id),ISM.Rel(id)];
            refListSpat.SE_GN1(kk,:)   = [refListDet(kk,1)-refWin(1),refListDet(kk,1)+refWin(2)];
        end

        if ISM.N(id) == 2
            refListSpat.GN2(kk) = id;
            refListSpat.AzEl_GN2(kk,:) = [ISM.Raz(id),ISM.Rel(id)];
            refListSpat.SE_GN2(kk,:)   = [refListDet(kk,1)-refWin(1),refListDet(kk,1)+refWin(2)];
        end    
    end
    
    %Finally select reflections for rendering
    selectIDs  = nan(length(refListSpat.toa),1);
    selectN    = nan(length(refListSpat.toa),1);
    selectAzEl = nan(length(refListSpat.toa),2);
    selectSE   = nan(length(refListSpat.toa),2);
    
    if psx.geo.prefN1 %Favor 1st order image sources!
    
        %Get 1st order ISM first
        N1Ren = ~isnan(refListSpat.N1); 

        selectIDs(N1Ren)    = refListSpat.N1(N1Ren);
        selectN(N1Ren)      = 1;
        selectAzEl(N1Ren,:) = refListSpat.AzEl_N1(N1Ren,:);
        selectSE(N1Ren,:)   = refListSpat.SE_N1(N1Ren,:);

        refToAssign = find(isnan(selectIDs));
        
        %Get 2nd order ISM if not already assigned
        if ~isempty(refToAssign)
            N2Ren = refListSpat.N2(refToAssign);

            selectIDs(refToAssign) = N2Ren;
            selectN(refToAssign)   = 2;
            selectAzEl(refToAssign,:) = refListSpat.AzEl_N2(refToAssign,:);
            selectSE(refToAssign,:)   = refListSpat.SE_N2(refToAssign,:);

            refToAssign = find(isnan(selectIDs));
        end

        %Get global N1/N2 reflection if still not asigned
        if ~isempty(refToAssign)
            GN1Ren = refListSpat.GN1(refToAssign);

            if ~isnan(GN1Ren)
                selectIDs(refToAssign) = GN1Ren;
                selectN(refToAssign)   = 3;
                selectAzEl(refToAssign,:) = refListSpat.AzEl_GN1(refToAssign,:);
                selectSE(refToAssign,:)   = refListSpat.SE_GN1(refToAssign,:);
            else
                GN2Ren = refListSpat.GN2(refToAssign);

                selectIDs(refToAssign) = GN2Ren;
                selectN(refToAssign)   = 4;
                selectAzEl(refToAssign,:) = refListSpat.AzEl_GN2(refToAssign,:);
                selectSE(refToAssign,:)   = refListSpat.SE_GN2(refToAssign,:);
            end
        end
        
    end
    
    if ~psx.geo.prefN1 %Just get 1st or 2nd order reflection from global list
        
        selectIDs(~isnan(refListSpat.GN2)) = refListSpat.GN2(~isnan(refListSpat.GN2));   
        selectIDs(~isnan(refListSpat.GN1)) = refListSpat.GN1(~isnan(refListSpat.GN1));   

        selectN(~isnan(refListSpat.GN2)) = 2;
        selectN(~isnan(refListSpat.GN1)) = 1;

        selectAzEl(~isnan(refListSpat.GN2),:) = refListSpat.AzEl_GN2(~isnan(refListSpat.GN2),:);
        selectAzEl(~isnan(refListSpat.GN1),:) = refListSpat.AzEl_GN1(~isnan(refListSpat.GN1),:);

        selectSE(~isnan(refListSpat.GN2),:)   = refListSpat.SE_GN2(~isnan(refListSpat.GN2),:);
        selectSE(~isnan(refListSpat.GN1),:)   = refListSpat.SE_GN1(~isnan(refListSpat.GN1),:);
    
    end
    
    
    %-------------------------------------%
    %Write final list with selected reflections
    refListSpat.selectIDs  = selectIDs;
    refListSpat.selectN    = selectN;
    refListSpat.selectAzEl = selectAzEl;
    refListSpat.selectAzEl_SH = selectAzEl;
    refListSpat.selectAzEl_SH(:,2) = 90-refListSpat.selectAzEl_SH(:,2);
    refListSpat.selectSE   = selectSE;
    refListSpat.selectSrcAzEl = [ISM.Saz(selectIDs),ISM.Sel(selectIDs)];
    refListSpat.selectSrcAzEl_SH = refListSpat.selectSrcAzEl;
    refListSpat.selectSrcAzEl_SH(:,2) = 90-refListSpat.selectSrcAzEl_SH(:,2);
    
    %Write refListSpat in psx struct
    psx.spat.refListSpat = refListSpat;
    
end

%% Assign pseudo-randomized DOAs if no ISM simulation 

%If no DOA information and no room dimensions available
if ~isfield(psx.geo,'doa') && ~psx.geo.roomDimensions
   
    %Add direct sound, could be from ISM too, but lets say we just measured
    %distance and azimuth/elevation in the room
    refListSpat.tDirect = tDirect;
    refListSpat.eDirect = eDirect;
    refListSpat.AzEl_Direct = sourceDirection;
    refListSpat.AzEl_Direct_SH = refListSpat.AzEl_Direct;
    refListSpat.AzEl_Direct_SH(:,2) = 90-refListSpat.AzEl_Direct_SH(:,2);
    refListSpat.SE_Direct = [tDirect-directWin(1),tDirect+directWin(2)];
    refListSpat.SrcAzEl_Direct = [sourceDirection(1), -sourceDirection(2)];%From the point of view of the source
    refListSpat.SrcAzEl_Direct_SH = refListSpat.SrcAzEl_Direct;
    refListSpat.SrcAzEl_Direct_SH(:,2) = 90 - refListSpat.SrcAzEl_Direct_SH(:,2);
    refListSpat.srcView = [mod(180+sourceDirection(:,1), 360) -sourceDirection(:,2)];
    
    %Save in struct
    refListSpat.toa = refListDet(1:nRef,1); %TOA/E from detection list
    refListSpat.e   = refListDet(1:nRef,2);
    
    %Just take pre-defined list and write in selectAzEl
    arbRefPattern = psx.geo.arbRefPattern;
    %Adjust according to source direction
    if sum(abs(sourceDirection)) ~= 0
        %[arbRefPattern(:,1),arbRefPattern(:,2)] = AKroomSimulationRotation(arbRefPattern(:,1),arbRefPattern(:,2),360-sourceDirection(1),sourceDirection(2));
        arbRefPattern(:,1) = mod(arbRefPattern(:,1) + sourceDirection(1),360);
        arbRefPattern(:,2) = arbRefPattern(:,2)-sourceDirection(2);
    end
    if size(arbRefPattern,1) < nRef
        error('Please define at least as many arbitrary reflection directions as intended to be spatialized according to variable nRef');
    end
    
    refListSpat.selectAzEl = arbRefPattern(1:nRef,:);
    refListSpat.selectAzEl_SH = arbRefPattern(1:nRef,:);
    refListSpat.selectAzEl_SH(:,2) = 90-refListSpat.selectAzEl_SH(:,2);
    refListSpat.selectSE   = [refListSpat.toa-refWin(1),refListSpat.toa+refWin(2)];

    %Write refListSpat in psx struct
    psx.spat.refListSpat = refListSpat;
    
end

%% Use predetermined DOAs, e.g., from SDM, if available

if isfield(psx.geo,'doa')
    
    %Convert doa vector to spherical coordinates
    [doa(:,1),doa(:,2),~] = cart2sph(psx.geo.doa(:,1),psx.geo.doa(:,2),psx.geo.doa(:,3));
    doa(:,1) = doa(:,1)/pi*180; doa(:,2) = doa(:,2)/pi*180;
    doa(:,1) = mod(doa(:,1),360);
    
    %Add direct sound information. tDirect is simply based on sourceDistance
    refListSpat.tDirect = tDirect;
    refListSpat.eDirect = eDirect;
    %Use value of DOA vector at this specific point in time. As the DOA
    %vector is already smoothed and windowed, result should be fine.
    %Averaging over a specific window could also be tested.
    sourceDirection = doa(tDirect,:);
    psx.geo.sourceDirection = sourceDirection; %Write in psx.geo as required later
    refListSpat.AzEl_Direct = sourceDirection;
    refListSpat.AzEl_Direct_SH = refListSpat.AzEl_Direct;
    refListSpat.AzEl_Direct_SH(:,2) = 90-refListSpat.AzEl_Direct_SH(:,2);
    refListSpat.SE_Direct = [tDirect-directWin(1),tDirect+directWin(2)];
    refListSpat.SrcAzEl_Direct = [sourceDirection(1), -sourceDirection(2)];%From the point of view of the source
    refListSpat.SrcAzEl_Direct_SH = refListSpat.SrcAzEl_Direct;
    refListSpat.SrcAzEl_Direct_SH(:,2) = 90 - refListSpat.SrcAzEl_Direct_SH(:,2);
    refListSpat.srcView = [mod(180+sourceDirection(:,1), 360) -sourceDirection(:,2)];
    
    %Save in struct
    refListSpat.toa = refListDet(1:nRef,1); %TOA/E from detection list
    refListSpat.e   = refListDet(1:nRef,2);
    
    %Go through TOA list and take respective values from DOA vector
    %Again, averad values from a small window could have advantages
    %In this test, the DOAs are stable arround the TOA values and thus
    %there would not be a big difference
    refListSpat.selectAzEl = doa(refListSpat.toa,:);
    refListSpat.selectAzEl_SH = doa(refListSpat.toa,:);
    refListSpat.selectAzEl_SH(:,2) = 90-refListSpat.selectAzEl_SH(:,2);
    refListSpat.selectSE   = [refListSpat.toa-refWin(1),refListSpat.toa+refWin(2)];

    %Write refListSpat in psx struct
    psx.spat.refListSpat = refListSpat;
    
end
end