clc, clearvars, close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs for NExT-ERA identification %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -----------------------------------------------------------------------
%%% Inputs
data = load('data_1.txt');
IsTimeVectorIncluded = true;
%%% Properties for SigPro -------------------------------------------------
fs    = 600; % (Hz) Sampling frequency
fr    = 200; % (Hz) Resampling frequency or fs
ffi   = 0.1; % (Hz) Cutoff freq. for high-pass filter
fff   = 50; % (Hz) Cutoff freq. for low-pass filter
Wndw  = []; % window indicated as point number, [] to inlcude all.
Trend = 1; % Apply detren (1) or not (0).
pAcel = 0; % Plot time-domain accel. (1) or not (0).
%%% Other properties ------------------------------------------------------
RefChannel = 1; % Reference channel indicated as scalar or vector.
ModelOrderRange  = 1:4:40;
NumClusters = 4;
MAC_Treshold = 0.95;
PLOT_MAC = 1;
PLOT_FFT = 1;
%%% Global channels -------------------------------------------------------
%%% En caso que se tengan los sensores _023 (cada sensor con las
%%% direcciones), entonces GlobalChannels = [1:3,7:12];
%%% 
%%% A continuación se proveen varias alternativas para indicar el vector.
%%%
%%% Indicar directamente o dejar vacío:
GlobalChannels = [];
%%% Indicar Sensores incluidos y número de canales por sensor (3):
% Sensores = [1,3];
% NumCanalesPorSensor = 3;
% GlobalChannels = ((Sensores(:)-1)*NumCanalesPorSensor*ones(1,NumCanalesPorSensor)+...
%     ((1:NumCanalesPorSensor)'*ones(1,length(Sensores)))')';
% GlobalChannels = GlobalChannels(:);
%%%
%%% Frecuencias proximadas ------------------------------------------------
FreqAprox = []; % (Hz) 
% FreqAprox = [...
%      4.40,  4.80;
%     17.50, 19.00;
%     39.00, 42.00;
%     69.00, 75.00;
%     ]; % (Hz)
%%% Call function ---------------------------------------------------------
%%% -----------------------------------------------------------------------
[P50,ID,GlobalChannels] = ID_with_NExT_ERA(...
    data,...
    IsTimeVectorIncluded,...
    fs,...
    fr,...
    ffi,...
    fff,...
    Wndw,...
    Trend,...
    pAcel,...
    RefChannel,...
    ModelOrderRange,NumClusters,MAC_Treshold,...
    PLOT_MAC,...
    PLOT_FFT,...
    FreqAprox,...
    GlobalChannels);
