% Copyright (c) Facebook, Inc. and its affiliates.
%
%% -- Paraspax --
%  -- Parameterization, Spatialization, and Extrapolation 
%     of Monaural Room Impulse Responses --
%
% function f = paraspaxPlot(x,plotMode,fs,figSize)
%
% Extensive plot function to plat data in x according to plotMode
%
% Output:
% -                     - 
%
% Input:        
% x                     - Data
% plotMode              - Plot mode
% fs                    - Sample Rate
%                         Default: 48kHz
% figSize               - [x,y] array with width/height of figure
%                         Default: []
%
% Dependencies: -
%   
% Code written 2019/2020 by JMA, Johannes M. Arend.


function f = paraspaxPlot(x,plotMode,fs,figSize)

if nargin < 3 || isempty(fs)
    fs =  48000;
end

if nargin < 4 || isempty(figSize)
    figSize = [600 340];
end

%% Define global plot variables
blue = [109/255, 132/255, 180/255];
orange = [234/255 151/255 114/255];
yellow = [214/255, 201/255, 114/255];
pink = [201/255, 122/255, 227/255];
green = [49/255, 162/255, 147/255];
greenLight = [142/255, 222/255, 205/255];
purple = [158/255, 119/255, 189/255];
red = [217/255, 70/255, 62/255];
gray = [150/255, 150/255, 150/255];
gray2 = [75/255, 75/255, 75/255];
brown = [160/255, 148/255, 132/255];

lineWidth = 1.5;
lineWidthRir = 1.2;
fontSize = 14;

lineWidthPaper = 2.5;
lineWidthRirPaper = 2.0;
lineWidthRirPaper2 = 1.5;
fontSizePaper = 20;
markerSizePaper = 8;

%% 
if strcmp(plotMode,'rir')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs, length(x));
    plot(tVec,x,'Linewidth',lineWidthRir,'Color',blue)
    
    xlim([0,tVec(end)]);
    ylim([min(min(x)),max(max(x))]);
    
    xlabel('Time [s]');
    ylabel('Amplitude');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%% 
if strcmp(plotMode,'rir2c')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs, length(x));
    plot(tVec,x(:,1),'Linewidth',lineWidthRir,'Color',blue)
    hold on;
    plot(tVec,x(:,2),'Linewidth',lineWidthRir,'Color',orange)
    
    xlim([0,tVec(end)]);
    ylim([min(min(x)),max(max(x))]);
    
    xlabel('Time [s]');
    ylabel('Amplitude');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%% 
if strcmp(plotMode,'rir2c_ms')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs*1000, length(x));
    plot(tVec,x(:,1),'Linewidth',lineWidthRir,'Color',blue)
    hold on;
    plot(tVec,x(:,2),'Linewidth',lineWidthRir,'Color',orange)
    
    xlim([0,tVec(end)]);
    ylim([min(min(x)),max(max(x))]);
    
    xlabel('Time in ms');
    ylabel('Amplitude');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%% 
if strcmp(plotMode,'rir2c_ms_paper')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs*1000, length(x));
    plot(tVec,x(:,1),'Linewidth',lineWidthRirPaper,'Color',blue)
    hold on;
    plot(tVec,x(:,2),'Linewidth',lineWidthRirPaper,'Color',orange)
    
    xlim([0,tVec(end)]);
    ylim([min(min(x)),max(max(x))]);
    
    xlabel('Time in ms');
    ylabel('Amplitude');
    
    set(gca,'FontSize',fontSizePaper);
    
    grid on;
    
end

%%
if strcmp(plotMode,'rirLog')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs, length(x));
    plot(tVec,20*log10(abs(x)),'Linewidth',lineWidthRir,'Color',blue)
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(x)))),max(max(20*log10(abs(x))))]);
   
    xlabel('Time [s]');
    ylabel('Amplitude [dB]');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'rirLog2c')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs, length(x));
    plot(tVec,20*log10(abs(x(:,1))),'Linewidth',lineWidthRir,'Color',blue)
    hold on;
    plot(tVec,20*log10(abs(x(:,2))),'Linewidth',lineWidthRir,'Color',orange)
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(x)))),max(max(20*log10(abs(x))))]);
   
    xlabel('Time [s]');
    ylabel('Amplitude [dB]');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'rirLog2c_paper')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(x)/fs*1000, length(x));
    plot(tVec,20*log10(abs(x(:,1))),'Linewidth',lineWidthRirPaper,'Color',blue)
    hold on;
    plot(tVec,20*log10(abs(x(:,2))),'Linewidth',lineWidthRirPaper,'Color',orange)
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(x)))),max(max(20*log10(abs(x))))]);
   
    xlabel('Time in ms');
    ylabel('Amplitude in dB');
    
    set(gca,'FontSize',fontSizePaper);
    
    grid on;
    
end


%%
if strcmp(plotMode,'rirDetRef')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    wf = x.refDet.wf;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    %tDirect = x.par.mon.rirOnset_20dB;
    fs = x.fs;
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    
    plot(tVec,abs(rir(1:mts2)),'Linewidth',lineWidthRir,'Color',blue)
    hold on;
    plot(tVec,wf,'Linewidth',lineWidth,'Color',yellow) %WF
    plot(tVec(refListDet(1:nRef,1)),abs(rir(refListDet(1:nRef,1))),'o','Linewidth',lineWidth,'Color',pink) %Peaks
    plot(tVec(refListDet(1:nRef,1)),abs(refListDet(1:nRef,2)),'s','Linewidth',lineWidth,'Color',orange) %Energy
    plot(ones(128,1)*tVec(mts2/2),linspace(0,1,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),0.88,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    legend('RIR','WF','Peak','RMS');
    
    xlim([0, tVec(end)])
    ylim([min(abs(rir)),max(abs(rir))]);
    
    xlabel('Time [ms]');
    ylabel('Absolute Amplitude');
    
    set(gca,'FontSize',fontSize);
    
end

%%
if strcmp(plotMode,'rirDetRef_paper')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    wf = x.refDet.wf;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    %tDirect = x.par.mon.rirOnset_20dB;
    fs = x.fs;
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    
    plot(tVec,abs(rir(1:mts2)),'Linewidth',lineWidthRirPaper,'Color',blue)
    hold on;
    plot(tVec,wf,'Linewidth',lineWidthPaper,'Color',yellow) %WF
    plot(tVec(refListDet(1:nRef,1)),abs(rir(refListDet(1:nRef,1))),'o','Linewidth',lineWidthPaper,'Color',pink,'Markersize',markerSizePaper) %Peaks
    plot(tVec(refListDet(1:nRef,1)),abs(refListDet(1:nRef,2)),'s','Linewidth',lineWidthPaper,'Color',orange,'Markersize',markerSizePaper) %Energy
    plot(ones(128,1)*tVec(mts2/2),linspace(0,1,128),'--','Linewidth',lineWidthPaper,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),0.915,['MT = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSizePaper,'Color',gray);
    
    legend('RIR','WF','Peak','RMS','Location','SouthWest');
    
    xlim([0, tVec(end)])
    ylim([min(abs(rir)),max(abs(rir))]);
    
    xlabel('Time in ms');
    ylabel('Amplitude');
    
    grid on;
    
    set(gca,'FontSize',fontSizePaper);
    
end

%% 
if strcmp(plotMode,'rirLogDetRef')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    %wf = x.refDet.wf;
    percTcurve = x.refDet.percTcurve;
    percTcurveLow = x.refDet.percTcurveLow;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    %tDirect = x.par.mon.rirOnset_20dB;
    fs = x.fs;
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    range = 60;
    
    p1 = plot(tVec,20*log10(abs(rir(1:mts2))),'Linewidth',lineWidthRir,'Color',blue);
    hold on;
    p2 = plot(tVec,percTcurve,'LineStyle','--','Linewidth',lineWidth,'Color',greenLight); %pThresh
    p3 = plot(tVec,percTcurveLow,'LineStyle','-','Linewidth',lineWidth,'Color',greenLight); %pThresh
    p4 = plot(tVec(refListDet(1:nRef,1)),20*log10(abs(rir(refListDet(1:nRef,1)))),'o','Linewidth',lineWidth,'Color',pink); %Peaks
    p5 = plot(tVec(refListDet(1:nRef,1)),20*log10(abs(refListDet(1:nRef,2))),'s','Linewidth',lineWidth,'Color',orange); %Energy
    plot(ones(128,1)*tVec(mts2/2),linspace(-range,0,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),max(20*log10(abs(rir)))-5,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    legend([p1,p2,p3,p4,p5],'RIR','pThresh','pThresh_l_o_w','Peak','RMS');
    
    xlim([0, tVec(end)])
    ylim([max(20*log10(abs(rir)))-range,max(20*log10(abs(rir)))]);
    
    xlabel('Time [ms]');
    ylabel('Magnitude [dB]');
    
    set(gca,'FontSize',fontSize);
    
end

%% 
if strcmp(plotMode,'rirLogDetRef_paper')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    %wf = x.refDet.wf;
    percTcurve = x.refDet.percTcurve;
    percTcurveLow = x.refDet.percTcurveLow;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    %tDirect = x.par.mon.rirOnset_20dB;
    fs = x.fs;
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    range = 60;
    
    p1 = plot(tVec,20*log10(abs(rir(1:mts2))),'Linewidth',lineWidthRirPaper,'Color',blue);
    hold on;
    p2 = plot(tVec,percTcurve,'LineStyle','--','Linewidth',lineWidthPaper,'Color',greenLight); %pThresh
    %p3 = plot(tVec,percTcurveLow,'LineStyle','-','Linewidth',lineWidth,'Color',greenLight); %pThresh
    p4 = plot(tVec(refListDet(1:nRef,1)),20*log10(abs(rir(refListDet(1:nRef,1)))),'o','Linewidth',lineWidthPaper,'Color',pink,'Markersize',markerSizePaper); %Peaks
    p5 = plot(tVec(refListDet(1:nRef,1)),20*log10(abs(refListDet(1:nRef,2))),'s','Linewidth',lineWidthPaper,'Color',orange,'MarkerSize',markerSizePaper); %Energy
    plot(ones(128,1)*tVec(mts2/2),linspace(-range,0,128),'--','Linewidth',lineWidthPaper,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),max(20*log10(abs(rir)))-5,['MT = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSizePaper,'Color',gray);
    
    legend([p1,p2,p4,p5],'RIR','mThresh','Peak','RMS','Location','SouthWest');
    
    xlim([0, tVec(end)])
    ylim([max(20*log10(abs(rir)))-range,max(20*log10(abs(rir)))]);
    
    xlabel('Time in ms');
    ylabel('Amplitude in dB');
    
    grid on;
    
    set(gca,'FontSize',fontSizePaper);
    
end

%% 
if strcmp(plotMode,'rirRT')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    semilogx(x.par.mon.freqVector,x.par.mon.TLundeby,'-*','Color',blue,'LineWidth',lineWidth,'MarkerSize',10);
    
    xlim([125,16000]);
    %ylim([min(x),max(x)]);
    
    xticks([125,250,500,1000,2000,4000,8000,16000])
    xticklabels({'125','250','500','1k','2k','4k','8k','16k'})
    
    legend('RT_6_0');
        
    xlabel('Frequency [Hz]');
    ylabel('Time [s]');
    
    grid on;
    set(gca,'FontSize',fontSize);
   
end

%%
if strcmp(plotMode,'rirMT')
    
    rir = x.rir;
    edpAbel = x.par.mon.edpAbel;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(rir)/fs, length(rir));
    plot(tVec,rir,'Linewidth',lineWidthRir,'Color',blue)
    hold on;
    plot(tVec,edpAbel,'Color',orange,'LineWidth',lineWidth);
    plot(ones(128,1)*tVec(mts2/2),linspace(min(rir),max(rir),128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),0.51,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    xlim([0,tVec(end)]);
    ylim([min(rir),max(edpAbel)]);
    
    legend('RIR','EDP_A_b_e_l','Location','east');
    
    xlabel('Time [s]');
    ylabel('Amplitude / Echo Density');
    
    set(gca,'FontSize',fontSize);
    
end

%% 
if strcmp(plotMode,'mag')

    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    fVec = linspace(0,fs/2,length(x));
    semilogx(fVec,20*log10(abs(x)),'Linewidth',lineWidth,'Color',blue)
    
    xlim([20,20000]);
    ylim([min(20*log10(abs(x))),max(20*log10(abs(x)))]);
    
    xticks([100,1000,10000])
    xticklabels({'100','1k','10k'})
    
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'refFilt')

    nRef = 3;
    
    refFilt = x.refDet.refFilt;
    fVec = x.refDet.fVecRefFilt;
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    semilogx(fVec,20*log10(abs(refFilt(:,1))),'Linewidth',lineWidth,'Color',blue);
    hold on;
    semilogx(fVec,20*log10(abs(refFilt(:,2))),'Linewidth',lineWidth,'Color',orange);
    semilogx(fVec,20*log10(abs(refFilt(:,3))),'Linewidth',lineWidth,'Color',pink);
    
    xlim([100,20000]);
    ylim([min(min(20*log10(abs(refFilt(:,1:nRef))))),max(max(20*log10(abs(refFilt(:,1:nRef)))))]);
    
    xticks([100,1000,10000])
    xticklabels({'100','1k','10k'})
    
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    
    legend('Reflection 1','Reflection 2','Reflection 3');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'refFilt_paper')

    nRef = 3;
    
    dsFilt = x.refDet.dsFilt;
    refFilt = x.refDet.refFilt;
    fVec = x.refDet.fVecRefFilt;
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    %3 Reflections
    p1 = semilogx(fVec,20*log10(abs(refFilt(:,1))),'Linewidth',lineWidthPaper,'Color',blue);
    hold on;
    p2 = semilogx(fVec,20*log10(abs(refFilt(:,2))),'Linewidth',lineWidthPaper,'Color',orange);
    p3 = semilogx(fVec,20*log10(abs(refFilt(:,3))),'Linewidth',lineWidthPaper,'Color',pink);
    
    xlim([100,20000]);
    ylim([min(min(20*log10(abs(refFilt(:,1:nRef))))),max(max(20*log10(abs(refFilt(:,1:nRef)))))]);
    
    xticks([100,1000,10000])
    xticklabels({'100','1k','10k'})
    
    xlabel('Frequency in Hz');
    ylabel('Magnitude in dB');
    
    legend('RF1','RF2','RF3');
    
    set(gca,'FontSize',fontSizePaper);
    
    grid on;
    
end

%% 
if strcmp(plotMode,'mag2c')

    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    fVec = linspace(0,fs/2,length(x));
    semilogx(fVec,20*log10(abs(x(:,1))),'Linewidth',lineWidth,'Color',blue)
    hold on;
    semilogx(fVec,20*log10(abs(x(:,2))),'Linewidth',lineWidth,'Color',orange)
    
    xlim([20,20000]);
    ylim([min(min(20*log10(abs(x)))),max(max(20*log10(abs(x))))]);
    
    xticks([100,1000,10000])
    xticklabels({'100','1k','10k'})
    
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%% 
if strcmp(plotMode,'ICfd')

    %Get rid of zeros as IC = 1 there
    onset = floor(min(AKonsetDetect(x,10,-20)));
    inWin = 32;
    if onset <= inWin
        x = [zeros(inWin*2,2);x];
        onset = onset+(inWin*2);
    end
    x = x(onset-inWin*2:end,:);
    x(:,1) = x(:,1).*supdeq_win(length(x),[inWin 0]);
    x(:,2) = x(:,2).*supdeq_win(length(x),[inWin 0]);
  
    blockSize = 512;
    nOct = 1;
    %Cut BRIR according to blocksize
    r = mod(length(x),blockSize);
    x = x(1:end-r,:);
    NFFT = length(x);
    hannWin = hann(blockSize);
    overLap = blockSize/2; %50 percent
    nBlocks = NFFT / (blockSize - overLap);
    IC_NFFT = 4096;

    icNum = zeros(nBlocks,IC_NFFT);
    icDenomHL = zeros(nBlocks,IC_NFFT);
    icDenomHR = zeros(nBlocks,IC_NFFT);
    for kk = 0:nBlocks-2

        ss = 1+(kk * overLap);
        se = ss+blockSize-1;

        brirPart = x(ss:se,:);
        brirPart(:,1) = brirPart(:,1) .* hannWin;
        brirPart(:,2) = brirPart(:,2) .* hannWin;

        HL = fft(brirPart(:,1),IC_NFFT);
        HR = fft(brirPart(:,2),IC_NFFT);

        icNum(kk+1,:) = HL .* conj(HR);
        icDenomHL(kk+1,:) = abs(HL).^2;
        icDenomHR(kk+1,:) = abs(HR).^2;

    end

    icNumFD = real(sum(icNum));
    icDenomFD = sqrt(sum(icDenomHL) .* sum(icDenomHR));
    icFD = icNumFD ./ icDenomFD;

    icFD = icFD(1:end/2+1);
    icFD = AKfractOctSmooth(icFD,'amp',fs,nOct);
    fVec = linspace(0,fs/2,IC_NFFT/2+1);

    %Actual plot
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    semilogx(fVec,icFD,'Linewidth',lineWidth,'Color',blue)
    
    xlim([20 20000]);
    %ylim([0 1]);
    
    xticks([100,1000,10000])
    xticklabels({'100','1k','10k'})
    
    xlabel('Frequency [Hz]');
    ylabel('Interaural Coherence');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'ICfd_paper')

    %Get rid of zeros as IC = 1 there
    onset = floor(min(AKonsetDetect(x,10,-20)));
    inWin = 32;
    if onset <= inWin
        x = [zeros(inWin*2,2);x];
        onset = onset+(inWin*2);
    end
    x = x(onset-inWin*2:end,:);
    x(:,1) = x(:,1).*supdeq_win(length(x),[inWin 0]);
    x(:,2) = x(:,2).*supdeq_win(length(x),[inWin 0]);
  
    blockSize = 512;
    nOct = 1;
    %Cut BRIR according to blocksize
    r = mod(length(x),blockSize);
    x = x(1:end-r,:);
    NFFT = length(x);
    hannWin = hann(blockSize);
    overLap = blockSize/2; %50 percent
    nBlocks = NFFT / (blockSize - overLap);
    IC_NFFT = 4096;

    icNum = zeros(nBlocks,IC_NFFT);
    icDenomHL = zeros(nBlocks,IC_NFFT);
    icDenomHR = zeros(nBlocks,IC_NFFT);
    for kk = 0:nBlocks-2

        ss = 1+(kk * overLap);
        se = ss+blockSize-1;

        brirPart = x(ss:se,:);
        brirPart(:,1) = brirPart(:,1) .* hannWin;
        brirPart(:,2) = brirPart(:,2) .* hannWin;

        HL = fft(brirPart(:,1),IC_NFFT);
        HR = fft(brirPart(:,2),IC_NFFT);

        icNum(kk+1,:) = HL .* conj(HR);
        icDenomHL(kk+1,:) = abs(HL).^2;
        icDenomHR(kk+1,:) = abs(HR).^2;

    end

    icNumFD = real(sum(icNum));
    icDenomFD = sqrt(sum(icDenomHL) .* sum(icDenomHR));
    icFD = icNumFD ./ icDenomFD;

    icFD = icFD(1:end/2+1);
    icFD = AKfractOctSmooth(icFD,'amp',fs,nOct);
    fVec = linspace(0,fs/2,IC_NFFT/2+1);

    %Actual plot
    %f = figure;
    %set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    semilogx(fVec,icFD,'Linewidth',lineWidthPaper,'Color',blue)
    
    xlim([20 20000]);
    %ylim([0 1]);
    
    xticks([100,1000,10000])
    xticklabels({'100','1k','10k'})
    
    xlabel('Frequency in Hz');
    ylabel('Interaural Coherence');
    title('Interaural Coherence');
    
    set(gca,'FontSize',fontSizePaper);
    
    grid on;
    
end

%%
if strcmp(plotMode,'ICfdtd')
    
    %Get rid of zeros as IC = 1 there
    onset = floor(min(AKonsetDetect(x,10,-20)));
    inWin = 32;
    x = x(onset-inWin*2:end,:);
    x(:,1) = x(:,1).*supdeq_win(length(x),[inWin 0]);
    x(:,2) = x(:,2).*supdeq_win(length(x),[inWin 0]);
    
    blockSize = 512;
    nOct = 3;
    %Cut BRIR according to blocksize
    r = mod(length(x),blockSize);
    x = x(1:end-r,:);
    NFFT = length(x);
    hannWin = hann(blockSize);
    nBlocks = NFFT / blockSize;
    IC_NFFT = 4096;

    startPoints = zeros(nBlocks-1,1);
    endPoints = zeros(nBlocks-1,1);
    icNum = zeros(nBlocks,IC_NFFT);
    icDenomHL = zeros(nBlocks,IC_NFFT);
    icDenomHR = zeros(nBlocks,IC_NFFT);
    for kk = 0:nBlocks-1

        ss = 1+(kk * blockSize);
        startPoints(kk+1) = ss;
        se = ss+blockSize-1;
        endPoints(kk+1) = se;

        brirPart = x(ss:se,:);
        brirPart(:,1) = brirPart(:,1) .* hannWin;
        brirPart(:,2) = brirPart(:,2) .* hannWin;

        HL = fft(brirPart(:,1),IC_NFFT);
        HR = fft(brirPart(:,2),IC_NFFT);

        icNum(kk+1,:) = HL .* conj(HR);
        icDenomHL(kk+1,:) = abs(HL).^2;
        icDenomHR(kk+1,:) = abs(HR).^2;

    end

    %Apply moving average with weighting
    l = -2:2;
    w = hann(length(l));
    zeroPad = zeros(abs(l(1)),IC_NFFT);
    icNumPad = [zeroPad;icNum;zeroPad];
    icDenomHLPad = [zeroPad;icDenomHL;zeroPad];
    icDenomHRPad = [zeroPad;icDenomHR;zeroPad];
    icNumW = zeros(size(icNumPad,1),size(icNumPad,2));
    icDenomHLW = zeros(size(icDenomHLPad,1),size(icDenomHLPad,2));
    icDenomHRW = zeros(size(icDenomHRPad,1),size(icDenomHRPad,2));
    for kk = abs(l(1)):nBlocks

        icNumW(kk+1,:) = sum(w.*icNumPad(kk+1+l(1):kk+1+l(end),:));
        icDenomHLW(kk+1,:) = sum(w.*icDenomHLPad(kk+1+l(1):kk+1+l(end),:));
        icDenomHRW(kk+1,:) = sum(w.*icDenomHRPad(kk+1+l(1):kk+1+l(end),:));

    end
    icNumW = icNumW(abs(l(1))+1:end-abs(l(1)),:);
    icDenomHLW = icDenomHLW(abs(l(1))+1:end-abs(l(1)),:);
    icDenomHRW = icDenomHRW(abs(l(1))+1:end-abs(l(1)),:);

    %icNumW = real(icNumW);
    %icDenomW = sqrt(icDenomHLW .* icDenomHRW);

    icTDFD = zeros(nBlocks,IC_NFFT);
    icTDFDVec = zeros(NFFT,IC_NFFT);
    for kk = 1:nBlocks-1
        icTemp = real(icNumW(kk,:) ./ sqrt(icDenomHLW(kk,:) .* icDenomHRW(kk,:)));
        if sum(isnan(icTemp)) == IC_NFFT
           icTemp = zeros(1,IC_NFFT); 
        end
        icTemp = AKfractOctSmooth(icTemp.','amp',fs,nOct); icTemp = icTemp.';
        icTDFD(kk,:) = icTemp;
        icTDFDVec(startPoints(kk):endPoints(kk),:) = icTDFDVec(startPoints(kk):endPoints(kk),:)+icTemp;
    end

    icTDFDVec = icTDFDVec(:,1:end/2+1);
    fVec = linspace(0,fs/2,IC_NFFT/2+1);
    
    %Actual plot
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    imagesc(1:size(icTDFDVec,1),fVec,icTDFDVec.');
    set(gca,'YDir','normal')
    ylim([0 20000]);
    xlim([0 NFFT]);

    xticksN = 0:fs/4:length(x);
    xticksN(1) = 1;
    xticklabelsA = xticksN/fs;
    xticklabelsA(1) = 0;
    set(gca,'Xtick',xticksN,'XTickLabel',xticklabelsA);
    yticksN = [100,1000,10000,20000];
    yticklabelsA = {'100','1k','10k','20k'};
    set(gca,'Ytick',yticksN,'YTickLabel',yticklabelsA);

    cb = colorbar;
    set(cb, 'ylim', [-1 1])
    set(cb, 'Ydir', 'normal')
    ylabel(cb, 'Interaural Coherence');
    colormap(flip(AKcolormaps('RdYlBu')));
    %colormap(jet);
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    
    set(gca,'FontSize',fontSize);
    
end

%%
if strcmp(plotMode,'monRefMap')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    tDirect = x.spat.refListSpat.tDirect;
    eDirect = x.spat.refListSpat.eDirect;
    pDirect = max(abs(x.rir));
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;

    s1 = scatter(tVec(tDirect),0,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    scatter(tVec(refListDet(:,1)),zeros(length(refListDet(:,1)),1),markerSize1,20*log10(abs(refListDet(:,3))),'filled');
    
    %Mark choosen reflections
    s2 = scatter(tVec(tDirect),0,markerSize3,'MarkerEdgeColor',pink);
    s3 = scatter(tVec(refListDet(1:nRef,1)),zeros(length(refListDet(1:nRef,1)),1),markerSize2,'MarkerEdgeColor',red);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-1,1,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),0.88,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    xlim([0, tVec(end)])
    ylim([-1 1]);
    
    yticks([])
    
    xlabel('Time [ms]');
    ylabel(cb, 'Amplitude [dB]')
    
    set(gca,'FontSize',fontSize);
    
    legend([s1,s2,s3],'Detected Reflections','Direct Sound','Selected Reflections');
    
    grid on;
    
end

%%
if strcmp(plotMode,'monRefMap_paper')
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    tDirect = x.spat.refListSpat.tDirect;
    eDirect = x.spat.refListSpat.eDirect;
    pDirect = max(abs(x.rir));
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;
    markerLineWidth = 2;

    s1 = scatter(tVec(tDirect),0,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    scatter(tVec(refListDet(:,1)),zeros(length(refListDet(:,1)),1),markerSize1,20*log10(abs(refListDet(:,3))),'filled');
    
    %Mark choosen reflections
    s2 = scatter(tVec(tDirect),0,markerSize3,'MarkerEdgeColor',pink,'LineWidth',markerLineWidth);
    s3 = scatter(tVec(refListDet(1:nRef,1)),zeros(length(refListDet(1:nRef,1)),1),markerSize2,'MarkerEdgeColor',red,'LineWidth',markerLineWidth);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-1,1,128),'--','Linewidth',lineWidthPaper,'Color',gray); %Mixing Time
    text(tVec(mts2/2+30),0.88,['MT = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSizePaper,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [-30 0]) %[round(20*log10(abs(min(refListDet(:,3)))))-5 0]
    %xlim([0, tVec(end)])
    xlim([0, 70])
    ylim([-1 1]);
    
    yticks([])
    
    xlabel('Time in ms');
    ylabel(cb, 'Amplitude in dB')
    
    set(gca,'FontSize',fontSizePaper);
    
    legend([s1,s2,s3],'Detected Reflections','Direct Sound','Selected Reflections','Location','NorthWest');
    
    grid on;
    
end

%%
if strcmp(plotMode,'spatRefMapAz')
    
    %Azimuth
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    tDirect = x.spat.refListSpat.tDirect;
    azDirect = x.spat.refListSpat.AzEl_Direct(:,1);
    elDirect = 90-x.spat.refListSpat.AzEl_Direct(:,2);
    if azDirect > 180
        azDirect = azDirect-360;
    end
    azRef = x.spat.refListSpat.selectAzEl(:,1);
    azRef(azRef>180) = azRef(azRef>180)-360;
    elRef = 90-x.spat.refListSpat.selectAzEl(:,2);
    pDirect = max(abs(x.rir));
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;

    s1 = scatter(tVec(tDirect),azDirect,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    s2 = scatter(tVec(refListDet(1:nRef,1)),azRef,markerSize1,20*log10(abs(refListDet(1:nRef,3))),'filled');
    %Dummy to generate minimum amplitude
    scatter(tVec(end)+10,0,markerSize1,20*log10(abs(refListDet(end,3))),'filled');
    
    %Mark choosen reflections
    s3 = scatter(tVec(tDirect),azDirect,markerSize3,'MarkerEdgeColor',pink);
    s4 = scatter(tVec(refListDet(1:nRef,1)),azRef,markerSize2,'MarkerEdgeColor',red);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-210,210,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),-150,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    xlim([0, tVec(end)])
    ylim([-180 180]);
    yticks([-180:60:180]);
    
    %Invert y axis
    set(gca, 'YDir','reverse')
    
    xlabel('Time [ms]');
    ylabel('Azimuth [°]');
    ylabel(cb, 'Amplitude [dB]')
    
    set(gca,'FontSize',fontSize);
    
    legend([s1,s3,s4],'Detected Reflections','Direct Sound','Selected Reflections');
    
    grid on;
    
end

%%
if strcmp(plotMode,'spatRefMapAz_paper')
    
    %Azimuth
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    tDirect = x.spat.refListSpat.tDirect;
    azDirect = x.spat.refListSpat.AzEl_Direct(:,1);
    elDirect = 90-x.spat.refListSpat.AzEl_Direct(:,2);
    if azDirect > 180
        azDirect = azDirect-360;
    end
    azRef = x.spat.refListSpat.selectAzEl(:,1);
    azRef(azRef>180) = azRef(azRef>180)-360;
    elRef = 90-x.spat.refListSpat.selectAzEl(:,2);
    pDirect = max(abs(x.rir));
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;
    markerLineWidth = 2;

    s1 = scatter(tVec(tDirect),azDirect,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    s2 = scatter(tVec(refListDet(1:nRef,1)),azRef,markerSize1,20*log10(abs(refListDet(1:nRef,3))),'filled');
    %Dummy to generate minimum amplitude
    scatter(tVec(end)+10,0,markerSize1,20*log10(abs(refListDet(end,3))),'filled');
    
    %Mark choosen reflections
    s3 = scatter(tVec(tDirect),azDirect,markerSize3,'MarkerEdgeColor',pink,'LineWidth',markerLineWidth);
    s4 = scatter(tVec(refListDet(1:nRef,1)),azRef,markerSize2,'MarkerEdgeColor',red,'LineWidth',markerLineWidth);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-210,210,128),'--','Linewidth',lineWidthPaper,'Color',gray); %Mixing Time
    text(tVec(mts2/2+30),-150,['MT = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSizePaper,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [-30 0]);%[round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    %xlim([0, tVec(end)])
    xlim([0, 70])
    ylim([-180 180]);
    yticks([-180:60:180]);
    
    %Invert y axis
    set(gca, 'YDir','reverse')
    
    xlabel('Time in ms');
    ylabel('Azimuth in degree');
    ylabel(cb, 'Amplitude in dB')
    
    set(gca,'FontSize',fontSizePaper);
    
    legend([s1,s3,s4],'Detected Reflections','Direct Sound','Selected Reflections','Location','NorthWest');
    
    grid on;
    
end
    
%%
if strcmp(plotMode,'spatRefMapEl')
    
    %Elevation
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    tDirect = x.spat.refListSpat.tDirect;
    azDirect = x.spat.refListSpat.AzEl_Direct(:,1);
    elDirect = x.spat.refListSpat.AzEl_Direct(:,2);
    if azDirect > 180
        azDirect = azDirect-360;
    end
    azRef = x.spat.refListSpat.selectAzEl(:,1);
    azRef(azRef>180) = azRef(azRef>180)-360;
    elRef = x.spat.refListSpat.selectAzEl(:,2);
    pDirect = max(abs(x.rir));
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;

    s1 = scatter(tVec(tDirect),elDirect,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    s2 = scatter(tVec(refListDet(1:nRef,1)),elRef,markerSize1,20*log10(abs(refListDet(1:nRef,3))),'filled');
    %Dummy to generate minimum amplitude
    scatter(tVec(end)+10,0,markerSize1,20*log10(abs(refListDet(end,3))),'filled');
    
    %Mark choosen reflections
    s3 = scatter(tVec(tDirect),elDirect,markerSize3,'MarkerEdgeColor',pink);
    s4 = scatter(tVec(refListDet(1:nRef,1)),elRef,markerSize2,'MarkerEdgeColor',red);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-90,90,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),-75,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    xlim([0, tVec(end)])
    ylim([-90 90]);
    yticks([-90:30:90]);
    
    
    xlabel('Time [ms]');
    ylabel('Elevation [°]');
    ylabel(cb, 'Amplitude [dB]')
    
    set(gca,'FontSize',fontSize);
    
    legend([s1,s3,s4],'Detected Reflections','Direct Sound','Selected Reflections');
    
    grid on;
    
end

%%
if strcmp(plotMode,'spatRefMapEl_paper')
    
    %Elevation
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    nRef = x.refDet.nRef;
    tDirect = x.spat.refListSpat.tDirect;
    azDirect = x.spat.refListSpat.AzEl_Direct(:,1);
    elDirect = x.spat.refListSpat.AzEl_Direct(:,2);
    if azDirect > 180
        azDirect = azDirect-360;
    end
    azRef = x.spat.refListSpat.selectAzEl(:,1);
    azRef(azRef>180) = azRef(azRef>180)-360;
    elRef = x.spat.refListSpat.selectAzEl(:,2);
    pDirect = max(abs(x.rir));
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;
    markerLineWidth = 2;

    s1 = scatter(tVec(tDirect),elDirect,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    s2 = scatter(tVec(refListDet(1:nRef,1)),elRef,markerSize1,20*log10(abs(refListDet(1:nRef,3))),'filled');
    %Dummy to generate minimum amplitude
    scatter(tVec(end)+10,0,markerSize1,20*log10(abs(refListDet(end,3))),'filled');
    
    %Mark choosen reflections
    s3 = scatter(tVec(tDirect),elDirect,markerSize3,'MarkerEdgeColor',pink,'LineWidth',markerLineWidth);
    s4 = scatter(tVec(refListDet(1:nRef,1)),elRef,markerSize2,'MarkerEdgeColor',red,'LineWidth',markerLineWidth);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-90,90,128),'--','Linewidth',lineWidthPaper,'Color',gray); %Mixing Time
    text(tVec(mts2/2+30),75,['MT = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSizePaper,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [-30 0]);%[round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    %xlim([0, tVec(end)])
    xlim([0, 70])
    ylim([-90 90]);
    yticks([-90:30:90]);
    
    
    xlabel('Time in ms');
    ylabel('Elevation in degree');
    ylabel(cb, 'Amplitude in dB')
    
    set(gca,'FontSize',fontSizePaper);
    
    legend([s1,s3,s4],'Detected Reflections','Direct Sound','Selected Reflections','Location','NorthWest');
    
    grid on;
    
end


%%
if strcmp(plotMode,'spatRefMapAz_extrap')
    
    %Azimuth
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    toaRef = x.extrap.refListSpat.toa;
    nRef = x.refDet.nRef;
    tDirect = x.extrap.refListSpat.tDirect;
    azDirect = x.extrap.refListSpat.AzEl_Direct(:,1);
    elDirect = x.extrap.refListSpat.AzEl_Direct(:,2);
    if azDirect > 180
        azDirect = azDirect-360;
    end
    azRef = x.extrap.refListSpat.selectAzEl(:,1);
    azRef(azRef>180) = azRef(azRef>180)-360;
    elRef = x.extrap.refListSpat.selectAzEl(:,2);
    pDirect = max(abs(x.rir))* x.extrap.refListSpat.ampFactorDirect;
    ampFactor = x.extrap.refListSpat.ampFactor;
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;

    s1 = scatter(tVec(tDirect),azDirect,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    s2 = scatter(tVec(toaRef),azRef,markerSize1,20*log10(abs(refListDet(1:nRef,3).*ampFactor)),'filled');
    %Dummy to generate minimum amplitude
    scatter(tVec(end)+10,0,markerSize1,20*log10(abs(refListDet(end,3))),'filled');
    
    %Mark choosen reflections
    s3 = scatter(tVec(tDirect),azDirect,markerSize3,'MarkerEdgeColor',pink);
    s4 = scatter(tVec(toaRef),azRef,markerSize2,'MarkerEdgeColor',red);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-210,210,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),-150,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    xlim([0, tVec(end)])
    ylim([-180 180]);
    yticks([-180:60:180]);
    
    
    xlabel('Time [ms]');
    ylabel('Azimuth [°]');
    ylabel(cb, 'Amplitude [dB]')
    
    set(gca,'FontSize',fontSize);
    
    legend([s1,s3,s4],'Detected Reflections','Direct Sound','Selected Reflections');
    
    grid on;
    
end

%%
if strcmp(plotMode,'spatRefMapEl_extrap')
    
    %Azimuth
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    rir = x.rir;
    mts2 = round(x.par.mon.mtAbel/1000*x.fs)*2;
    refListDet = x.refDet.refListDet;
    toaRef = x.extrap.refListSpat.toa;
    nRef = x.refDet.nRef;
    tDirect = x.extrap.refListSpat.tDirect;
    azDirect = x.extrap.refListSpat.AzEl_Direct(:,1);
    elDirect = x.extrap.refListSpat.AzEl_Direct(:,2);
    if azDirect > 180
        azDirect = azDirect-360;
    end
    azRef = x.extrap.refListSpat.selectAzEl(:,1);
    azRef(azRef>180) = azRef(azRef>180)-360;
    elRef = x.extrap.refListSpat.selectAzEl(:,2);
    pDirect = max(abs(x.rir))* x.extrap.refListSpat.ampFactorDirect;
    ampFactor = x.extrap.refListSpat.ampFactor;
    fs = x.fs;
    cm = AKcolormaps('Blues');
    tVec = linspace(0, length(rir(1:mts2))/fs*1000, length(rir(1:mts2)));
    markerSize1 = 70;
    markerSize2 = 120;
    markerSize3 = 120;
    
    s1 = scatter(tVec(tDirect),elDirect,markerSize1,20*log10(abs(pDirect)),'filled');
    hold on;
    s2 = scatter(tVec(toaRef),elRef,markerSize1,20*log10(abs(refListDet(1:nRef,3).*ampFactor)),'filled');
    %Dummy to generate minimum amplitude
    scatter(tVec(end)+10,0,markerSize1,20*log10(abs(refListDet(end,3))),'filled');
    
    %Mark choosen reflections
    s3 = scatter(tVec(tDirect),elDirect,markerSize3,'MarkerEdgeColor',pink);
    s4 = scatter(tVec(toaRef),elRef,markerSize2,'MarkerEdgeColor',red);
    
    plot(ones(128,1)*tVec(mts2/2),linspace(-90,90,128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts2/2+50),-75,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    cb = colorbar;
    colormap(cm);
    
    set(cb, 'ylim', [round(20*log10(abs(min(refListDet(:,3)))))-5 0])
    xlim([0, tVec(end)])
    ylim([-90 90]);
    yticks([-90:30:90]);
    
    
    xlabel('Time [ms]');
    ylabel('Elevation [°]');
    ylabel(cb, 'Amplitude [dB]')
    
    set(gca,'FontSize',fontSize);
    
    legend([s1,s3,s4],'Detected Reflections','Direct Sound','Selected Reflections');
    
    grid on;
    
end

%% 
if strcmp(plotMode,'revLevel_EDC')
    
    rir = x.rir;
    mts = round(x.par.mon.mtAbel/1000*x.fs);
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(rir)/fs, length(rir));
   
    p1 = plot(tVec,20*log10(abs(rir)),'Linewidth',lineWidthRir,'Color',blue);
    hold on;
    p3 = plot(tVec(mts),20*log10(x.rev.edcMethod.ampAtMT),'s','Color',orange,'MarkerFaceColor',orange);
    p4 = plot(tVec,20*log10(x.rev.edcMethod.straightLine),'Linewidth',2,'Linestyle','-','Color',brown);
    p5 = plot(tVec,20*log10(x.rev.edcMethod.logLine),'Linewidth',2,'Linestyle','-','Color',purple);
    p6 = plot(tVec,20*log10(x.rev.edcMethod.linLine),'Linewidth',2,'Linestyle','-','Color',green);
    p7 = plot(tVec,20*log10(x.rev.edcMethod.interpLine),'Linewidth',2,'Linestyle','-','Color',orange);
    p2 = plot(tVec,20*log10(x.rev.edcMethod.edcLevel),'Linewidth',2,'Color',red);
    
    plot(ones(128,1)*tVec(mts),linspace(-180,max(max(20*log10(abs(rir)))),128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts+50),-4.5,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(rir)))),max(max(20*log10(abs(rir))))]);
   
    xlabel('Time [s]');
    ylabel('Amplitude [dB]');
    
    legend([p2,p7,p4,p5,p6],'Adjusted EDC','EDC Polyfit','MT Level','Logarithmic Fade','Linear Fade');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%% 
if strcmp(plotMode,'revLevel_RMS')
    
    rir = x.rir;
    mts = round(x.par.mon.mtAbel/1000*x.fs);
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(rir)/fs, length(rir));
   
    p1 = plot(tVec,20*log10(abs(rir)),'Linewidth',lineWidthRir,'Color',blue);
    hold on;
    p2 = plot(tVec,20*log10(x.rev.rmsMethod.rmsCurve),'Linewidth',lineWidthRir,'Color',red);
    p3 = plot(tVec(mts),20*log10(x.rev.rmsMethod.ampAtMT),'s','Color',orange,'MarkerFaceColor',orange);
    p4 = plot(tVec,20*log10(x.rev.rmsMethod.straightLine),'Linewidth',2,'Linestyle','-','Color',brown);
    p5 = plot(tVec,20*log10(x.rev.rmsMethod.logLine),'Linewidth',2,'Linestyle','-','Color',purple);
    p6 = plot(tVec,20*log10(x.rev.rmsMethod.linLine),'Linewidth',2,'Linestyle','-','Color',green);
    p7 = plot(tVec,20*log10(x.rev.rmsMethod.interpLine),'Linewidth',2,'Linestyle','-','Color',orange);
    
    plot(ones(128,1)*tVec(mts),linspace(-180,max(max(20*log10(abs(rir)))),128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts+50),-4.5,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(rir)))),max(max(20*log10(abs(rir))))]);
   
    xlabel('Time [s]');
    ylabel('Amplitude [dB]');
    
    legend([p2,p7,p4,p5,p6],'Sliding RMS','RMS Polyfit','MT Level','Logarithmic Fade','Linear Fade');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'revLevel_MAX')
    
    rir = x.rir;
    mts = round(x.par.mon.mtAbel/1000*x.fs);
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(rir)/fs, length(rir));
   
    p1 = plot(tVec,20*log10(abs(rir)),'Linewidth',lineWidthRir,'Color',blue);
    hold on;
    p2 = plot(tVec,20*log10(x.rev.maxMethod.maxCurve),'Linewidth',lineWidthRir,'Color',red);
    p3 = plot(tVec(mts),20*log10(x.rev.maxMethod.ampAtMT),'s','Color',orange,'MarkerFaceColor',orange);
    p4 = plot(tVec,20*log10(x.rev.maxMethod.straightLine),'Linewidth',2,'Linestyle','-','Color',brown);
    p5 = plot(tVec,20*log10(x.rev.maxMethod.logLine),'Linewidth',2,'Linestyle','-','Color',purple);
    p6 = plot(tVec,20*log10(x.rev.maxMethod.linLine),'Linewidth',2,'Linestyle','-','Color',green);
    p7 = plot(tVec,20*log10(x.rev.maxMethod.interpLine),'Linewidth',2,'Linestyle','-','Color',orange);
    
    plot(ones(128,1)*tVec(mts),linspace(-180,max(max(20*log10(abs(rir)))),128),'--','Linewidth',lineWidth,'Color',gray); %Mixing Time
    text(tVec(mts+50),-4.5,['MT_A_b_e_l = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSize,'Color',gray);
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(rir)))),max(max(20*log10(abs(rir))))]);
   
    xlabel('Time [s]');
    ylabel('Amplitude [dB]');
    
    legend([p2,p7,p4,p5,p6],'Sliding MAX','MAX Polyfit','MT Level','Logarithmic Fade','Linear Fade');
    
    set(gca,'FontSize',fontSize);
    
    grid on;
    
end

%%
if strcmp(plotMode,'revLevel_MAX_paper')
    
    rir = x.rir;
    mts = round(x.par.mon.mtAbel/1000*x.fs);
    refListDet = x.refDet.refListDet;
    
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    tVec = linspace(0, length(rir)/fs*1000, length(rir));
   
    p1 = plot(tVec,20*log10(abs(rir)),'Linewidth',lineWidthRirPaper2,'Color',blue);
    hold on;
    p2 = plot(tVec,20*log10(x.rev.maxMethod.maxCurve),'Linewidth',2,'Color',orange);
    p3 = plot(tVec(mts),20*log10(x.rev.maxMethod.ampAtMT),'s','Color',red,'MarkerFaceColor',red);
    %p4 = plot(tVec,20*log10(x.rev.maxMethod.straightLine),'Linewidth',2,'Linestyle','-','Color',brown);
    %p5 = plot(tVec,20*log10(x.rev.maxMethod.logLine),'Linewidth',2,'Linestyle','-','Color',purple);
    %p6 = plot(tVec,20*log10(x.rev.maxMethod.linLine),'Linewidth',2,'Linestyle','-','Color',green);
    p7 = plot(tVec,20*log10(x.rev.maxMethod.interpLine),'Linewidth',2.5,'Linestyle','-','Color',red);
    
    %Intersection point
    p8 = plot(tVec(refListDet(2,1)),20*log10(x.rev.maxMethod.interpLine(refListDet(2,1))),'x','Linewidth',lineWidthPaper,'Color',gray2,'MarkerSize',18); %Energy
    
    plot(ones(128,1)*tVec(mts),linspace(-180,max(max(20*log10(abs(rir)))),128),'--','Linewidth',lineWidthPaper,'Color',gray); %Mixing Time
    text(tVec(mts+50),-4.5,['MT = ',num2str(round(x.par.mon.mtAbel)),' ms'],'FontSize',fontSizePaper,'Color',gray);
    
    xlim([0,tVec(end)]);
    ylim([min(min(20*log10(abs(rir)))),max(max(20*log10(abs(rir))))]);
   
    xlabel('Time in ms');
    ylabel('Amplitude in dB');
    
    legend([p1,p2,p7,p8],'RIR','Sliding MAX','MAX Polyfit','Reverb Level');
    
    set(gca,'FontSize',fontSizePaper);
    
    grid on;
    
end

%%
if strcmp(plotMode,'sourceDirectivity')
    
    sourceDirectivity = x.brir.sourceDirectivity;
    filterLength = 4096;

    %Choose frequency band
    frac        = 1;   %1 - octave filter bank, 3 - third octave filter bank
    fc          = [250,1000,4000,8000]; %Center frequency
    filterOrder = 4;

    %Get corner frequencies
    if frac == 1 
        %fcentreOctave = 10^3 * (2 .^ [-6:4])
        q = 2^(1/2);
        for k = 1:length(fc)
            fupper(k) = fc(k) * q;
            flower(k) = fc(k) / q;
        end
    end

    if frac == 3
        %fcentreThirdOctave  = 10^3 * (2 .^ ([-18:13]/3))
        q = 2^(1/6);
        for k = 1:length(fc)
            fupper(k) = fc(k) * q;
            flower(k) = fc(k) / q;
        end
    end

    coordinates = getEquiangleGrid(1,90);
    coordinates(end+1,:) = [360,90]; %For polar plot
    directivityData = nan(size(coordinates,1),length(fc));

    for b = 1:length(fc)
        for k = 1:size(coordinates,1)

            lsDirectivity = AKisht(sourceDirectivity.SH.coeff, sourceDirectivity.SH.doFFT, [coordinates(k,1) coordinates(k,2)], sourceDirectivity.SH.SHTmode, sourceDirectivity.SH.isEven, sourceDirectivity.SH.compact, sourceDirectivity.SH.SHmode);
            lsDirectivity = real(lsDirectivity);
            lsDirectivity = AKinterpolateSpectrum(lsDirectivity, sourceDirectivity.SH.f, filterLength, {'nearest' 'linear' 'nearest'}, sourceDirectivity.SH.fs);

            %Find nan if there is nan and replace with nearest neighbour
            %Little tricky code from mathworks page, but it works and is shorter
            %than running through the array...
            m = flipud(lsDirectivity);
            t = ~isnan(m);
            ii = cumsum(t);
            ii(ii == 0) = 1;
            ii = bsxfun(@plus,[0,ii(end,1:end-1)],ii);
            m1 = m(t);
            lsDirectivity = flipud(m1(ii));
            clear m t ii m1;

            %Directivity Filter
            df = lsDirectivity;
            df = AKsingle2bothSidedSpectrum(df,true);
            df = real(ifft(df));
            df = circshift(df,filterLength/2-1);
            df = df .* hann(length(df));

            %Filter df 
            df = AKfilter(df,'bp',[flower(b) fupper(b)],0,x.fs,filterOrder,'butter');
  
            %Write in directivityData array
            directivityData(k,b) = 20*log10(sqrt(sum(df.^2)));

        end
    end

    %Normalize
    for b = 1:length(fc)
        directivityData(:,b) = directivityData(:,b)-max(directivityData(:,b));
    end

    %Plot
    f = figure;
    set(f,'units','points','position',[200,200,figSize(1),figSize(2)])
    
    dbRange = [-30 0];
    AKp(directivityData(:,1),'x6','az',coordinates(:,1),'dr', dbRange,'lw',2.5,'c',blue);
    hold on;
    AKp(directivityData(:,2),'x6','az',coordinates(:,1),'dr', dbRange,'lw',2.5,'c',orange);
    AKp(directivityData(:,3),'x6','az',coordinates(:,1),'dr', dbRange,'lw',2.5,'c',green);
    AKp(directivityData(:,4),'x6','az',coordinates(:,1),'dr', dbRange,'lw',2.5,'c',pink);

    set(gca,'FontSize',fontSize);
    set(gca,'LineWidth',1.5);
    legend('250 Hz','1 kHz','4 kHz','8 kHz');
    
end
end