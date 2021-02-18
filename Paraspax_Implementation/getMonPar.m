% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function par_omni = getMonPar(rir,fs,bandsPerOct)
%
% Function to calculate basic room acoustic parameters based on a single
% channel RIR (mostly omnidirectional). Most parameters are calculated
% using the ITA toolbox.
%
% Output:
% par_omni              - Struct with various room acoustic parameters
%
% Input:        
% rir                   - Single channel RIR
% fs                    - Sample rate of RIR
% bandsPerOct           - Bands per octave for analysis.
%                         Default: 3 (1/3 octave bands)
%
% Dependencies: ITA toolbox, AKtools
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function parOmni = getMonPar(rir,fs,bandsPerOct)

if size(rir,2) > size(rir,1)
    rir = rir.';
end

if nargin < 2
    error('Please specify sample rate of RIR');
end

if nargin < 3
    bandsPerOct = 3;
end

%% Get basic parameters with ITA toolbox

rirITA = var2ITA(rir,fs);

freqRange = [20 20000];
raResults = ita_roomacoustics(rirITA, 'EDT', 'T20', 'T30', 'T60', 'C50', 'C80', 'D50', 'D80', 'Center_Time', 'EDC', 'T_Lundeby', 'freqRange', freqRange, 'bandsPerOctave', bandsPerOct);
raResultsBB = ita_roomacoustics(rirITA, 'EDT', 'T20', 'T30', 'T60', 'C50', 'C80', 'D50', 'D80', 'Center_Time', 'EDC', 'T_Lundeby', 'freqRange', freqRange, 'broadbandAnalysis');

%% Calculate DRR

DRR = getDRR(rir,1,fs); %Paraspax according to window later used for energy detection in direct sound. Results in 1.5 ms window
DRRZ = getDRR(rir,2.5,fs); %Zahorik2002, 2.5 ms

%% Get mixing time

%According to Abel et al.
[mtAbel, edpAbel] = AKmixingTimeAbel(rir,1024,fs,0);

%Perceptual mixing time according to Lindau et al., based on Abel
[tmp50_data_based, tmp95_data_based] = AKdataBasedMixingTime(rir,1024,fs,0,false);

%% Save in struct

parOmni.freqRange = freqRange;
parOmni.freqVector = raResults.T20.freqVector;
parOmni.bandsPerOct = bandsPerOct;
parOmni.EDT = raResults.EDT.freqData;
parOmni.T20 = raResults.T20.freqData;
parOmni.T30 = raResults.T30.freqData;
parOmni.T60 = raResults.T60.freqData;
parOmni.TLundeby = raResults.T_Lundeby.freqData;
parOmni.C50 = raResults.C50.freqData;
parOmni.C80 = raResults.C80.freqData;
parOmni.D50 = raResults.D50.freqData;
parOmni.D80 = raResults.D80.freqData;
parOmni.centerTime = raResults.Center_Time.freqData;
parOmni.EDC = raResults.EDC.timeData; %Has to be shifted to fit the RIR (length(psx.par.mon.EDC_BB) + psx.par.mon.rirOnset_20dB)

%Various results
parOmni.DRR = DRR;
parOmni.DDR_Zahorik = DRRZ;
parOmni.mtAbel = mtAbel;
parOmni.edpAbel = edpAbel;
parOmni.mtLindau50 = tmp50_data_based;
parOmni.mtLindau95 = tmp95_data_based;
parOmni.rirOnset_20dB = floor(AKonsetDetect(rir,10,-20,'rel')); %Set to -20dB in accordance with ITA toolbox and Zahorik2002
parOmni.fs = fs;

%Broadband results
parOmni.EDT_BB = raResultsBB.EDT.freqData;
parOmni.T20_BB = raResultsBB.T20.freqData;
parOmni.T30_BB = raResultsBB.T30.freqData;
parOmni.T60_BB = raResultsBB.T60.freqData;
parOmni.TLundeby_BB = raResultsBB.T_Lundeby.freqData;
parOmni.C50_BB = raResultsBB.C50.freqData;
parOmni.C80_BB = raResultsBB.C80.freqData;
parOmni.D50_BB = raResultsBB.D50.freqData;
parOmni.D80_BB = raResultsBB.D80.freqData;
parOmni.centerTime_BB = raResultsBB.Center_Time.freqData;
parOmni.EDC_BB = raResultsBB.EDC.timeData; %Has to be shifted to fit the RIR (length(psx.par.mon.EDC_BB) + psx.par.mon.rirOnset_20dB). Or generate with 'method' 'noCut'...
parOmni.EDC_BB_AK = AKedc(rir,true); %EDC with AKtools (in dB). No time shift applied here. Is just around 0 till direct sound kicks in

end