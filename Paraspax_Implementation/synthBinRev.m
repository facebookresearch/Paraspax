% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function binRev = synthBinRev(rir,binNoise,blockSizeRir,blockSizeNoise)
%
% Function to synthesize diffuse binaural reverberation based on a monaural RIR. 
% The function applies segmentwise convolution of the monaural RIR with binaural
% white noise. For more details, see [1],[2].
%
% Output:
% binRev                - Binaural reverberation
%
% Input:        
% rir                   - Room impulse response
% binNoise              - Binaural white noise
% blockSizeRir          - Blocksize of rir segments used for convolution
% bockSizeNoise         - Bocksize of noise segments used for convolution
%
% Dependencies: -
%
% References:
% [1] C. Pörschmann, P. Stade, and J. M. Arend, Binauralization of omnidirectional 
% room impulse responses-algorithm and technical evaluation,
% Proc. 20th Int. Conf. Digit. Audio Eff. (DAFx-17), Edinburgh, UK, pp. 345-352, 2017.
%
% [2] J. M. Arend, S. V. Amengual Garí, C. Schissler, F. Klein, and P. W. Robinson, 
% “Six-Degrees-of-Freedom Parametric Spatial Audio Based on One Monaural Room Impulse Response,” 
% J. Audio Eng. Soc., vol. 69, no. 7/8, pp. 557–575, 2021. ﻿
% https://doi.org/10.17743/jaes.2021.0009
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function binRev = synthBinRev(rir,binNoise,blockSizeRir,blockSizeNoise)

%Init variables
binRev = zeros(length(rir) + blockSizeNoise,2);
startNoiseBlock = 1;
endNoiseBlock = blockSizeNoise;
winNoise = hann(blockSizeNoise);
winRir = hann(blockSizeRir);

for kk=1:blockSizeRir:length(rir)-(blockSizeRir-1)  
    
    %Update start and end points of noise block
    startNoiseBlock = startNoiseBlock + blockSizeNoise;
    endNoiseBlock   = endNoiseBlock + blockSizeNoise;

    % Randomize if noise is not long enough. Randomization so not the same
    % correlated blocks are applied twice for the convolution
    if endNoiseBlock>length(binNoise)
        startNoiseBlock=floor(rand(1)*blockSizeNoise)+1;
        endNoiseBlock=startNoiseBlock+blockSizeNoise-1;
    end
    
    %Get blocks for convolution in big loop, apply windows, and convolve
    %with monaural RIR, and perform overlapp add. All in once, because its fun...
    rirBlock = rir(kk:kk+blockSizeRir-1).*winRir;
    noiseBlock = binNoise(startNoiseBlock:endNoiseBlock,:).*winNoise;
    binRev(kk:kk+blockSizeRir+blockSizeNoise-2,1)=binRev(kk:kk+blockSizeRir+blockSizeNoise-2,1)+convFFT(rirBlock,noiseBlock(:,1));
    binRev(kk:kk+blockSizeRir+blockSizeNoise-2,2)=binRev(kk:kk+blockSizeRir+blockSizeNoise-2,2)+convFFT(rirBlock,noiseBlock(:,2));
end

%Cut binaural reverb according to lenght of monaural RIR
binRev=binRev(length(binRev)-length(rir)+1:length(binRev),:);

%Perform energy normalization based on energy in the middle of the RIRs
startNormBlock = floor(1/4*length(rir));
endNormBlock=floor(3/4*length(rir));
%Energy of rir
M1=sqrt(sum(rir(startNormBlock:endNormBlock).^2));
%Energy of diffuse BRIR for left and right channel
B1=sqrt(sum(binRev(startNormBlock:endNormBlock,:).^2));
%Apply normalization value with mean for left/right energy of binaural
%reverberation
binRev = binRev*(M1/mean(B1));

end