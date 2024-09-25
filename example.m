% Example of mixing and filtering white noise such that the output signals
% exhibit a specific spatial coherence. Optionally, it is possible to
% induce spectral smoothness, balance, or both jointly.
%
% The script works for any arbitrary 3-D microphone constellation.
% It is possible to mix different audio signals (replace the white noise).
% However, the input signals must exhibit equal power to generate
% an accurate spatial coherence.
%
% Dependencies
%       mix_signals.m
%       ccoherence.m
%
% Related paper
%       D. Mirabilii, S. J. Schlecht, E.A.P. Habets,
%       Generating coherence-constrained multisensor signals using
%       balanced mixing and spectrally smooth filters, The Journal
%       of the Acoustical Society of America, Vol. 149, 1425, 2021.
%
% Authors
%       Daniele Mirabilii
%       International Audio Laboratories of Erlangen, Germany
%       daniele.mirabilii@audiolabs-erlangen.de
%
%       Sebastian Jiro Schlecht
%       Aalto University, Finland
%       sebastian.schlecht@aalto.fi
%
%       Emanuël A.P. Habets
%       International Audio Laboratories of Erlangen, Germany
%       emanuel.habets@audiolabs-erlangen.de
%
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany
% Copyright (c) 2020 Sebastian J. Schlecht
%
% ----------------------------------------------------------------------------------
%
% This code is based on the existing code: https://github.com/ehabets/ANF-Generator
%
% E.A.P. Habets, I. Cohen and S. Gannot, 'Generating nonstationary
% multisensor signals under a spatial coherence constraint,'
% Journal of the Acoustical Society of America, Vol. 124, Issue 5,
% pp. 2911-2917, Nov. 2008.

close all
clear variables
clc

addpath('./functions');

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontSize',14)

% Initialization
Fs = 16000;                     % Sample frequency (Hz)
params.Fs = Fs;
K = 1024;                       % FFT length
params.K = K;
sc_type = 'spherical';          % Noise-field coherence model: 'corcos', 'spherical', 'cylindrical'
decomposition = 'EVD';          % Type of decomposition: 'EVD' or 'CHD'
processing = 'balanced+smooth'; % Processing method: 'standard', 'smooth', 'balanced', 'balanced+smooth'
dur = 10;                       % Input duration in seconds
L = dur*Fs;                     % Data length

% Additional parameter for the Corcos model
params.speed = 20;              % km/h
params.direction = 60;          % Degree w.r.t. "North" (y-axis) [anti-clockwise]

% Sensors position (arbitrary 2/3-D array) xyz in [m]
m1 = [0.12, 0, 0];              % First sensor coordinates
m2 = [0.08, 0, 0];              % Second sensor coordinates
m3 = [0.04, 0, 0];              % Third sensor coordinates
m4 = [0, 0, 0];                 % Fourth sensor coordinates
mm = [m1;m2;m3;m4];
M = length(mm(:,1));            % Number of channels
params.mm = mm;

% Summary of parameters
fprintf('Number of channels: %d\n',M)
fprintf('Spatial coherence: %s\n',sc_type)
fprintf('Decomposition: %s\n',decomposition)
fprintf('Processing: %s\n\n',processing)

% Generate target spatial coherence
DC = generate_target_coherence(sc_type,params);

% Generate mixing matrix with target spatial coherence (optional balanced/smooth)
[C, C_none] = mixing_matrix(DC,decomposition,processing);

% Evaluate balance and smoothness before and after applying the chosen processing method
evaluate_mixing_matrix(C_none,C,decomposition,processing,params,sc_type);

% Generate M mutually independent input signals of length L (replace with real audio at will)
n = randn(L,M);

% Generate sensor signals with target spatial coherence
x = mix_signals(n,C);

% Compare target and generated coherence

% Estimate generate coherence from the output signals
DC_gen = mccoherence(x,K,K/4);

% Compute coherence error (between target and generated signal coherence)
xi = sum(sum(abs(DC_gen - DC).^2,1),2);
xi_avg = pow2db(mean(xi));

% Plot example of generated vs target spatial coherence
Freqs = linspace(0,Fs/2,K/2+1); %vector of frequencies in Hz
switch sc_type
    case {'spherical','cylindrical'}
        % Case diffuse (real-valued coherence)
        fig = figure();
        for m =1:M-1
            subplot(M-1,1,m);
            plot(Freqs/1000,real(squeeze(DC_gen(1,m+1,:))),'-k','LineWidth',2)
            hold on;
            plot(Freqs/1000,real(squeeze(DC(1,m+1,:))),'-.b','LineWidth',2)
            hold off;
            title(sprintf('Sensors 1-%d',m+1));
            axis([0 Fs/2000 -1 1]);
            if m==1
                lgd = legend(sprintf('%s %s:\n \\xi = %2.1f dB',decomposition,processing,xi_avg),'Target');
                lgd.Location = 'northeast';
            end
            if m~=(M-1)
                set(gca,'XTickLabel',[]);
            end
            grid on;
        end
        han=axes(fig,'visible','off');
        han.YLabel.Visible='on';
        han.YLabel.Position(1) = -0.1;
        han.XLabel.Position(2) = -0.05;
        han.XLabel.Visible='on';
        ylabel(han,'Spatial coherence');
        xlabel(han,'frequency [kHz]');

    otherwise
        % Case 'corcos' or generally complex-valued coherence
        fig = figure();
        for m =1:(M-1)
            subp1 = subplot(2*(M-1),1,(2*m -1));
            plot(Freqs/1000,real(squeeze(DC_gen(1,m+1,:))),'-k','LineWidth',2);
            hold on;
            plot(Freqs/1000,real(squeeze(DC(1,m+1,:))),'-.b','LineWidth',2)
            hold off;
            grid on;
            subp1.Position(2) = subp1.Position(2);
            axis([0 Fs/2000  -1 1]);
            ylabel('Real')
            title(sprintf('Sensors 1-%d',m+1));
            if m==1
                lgd = legend(sprintf('%s %s:\n \\xi = %2.1f dB',decomposition,processing,xi_avg),'Target');
                lgd.Location = 'northeast';
            end
            set(gca,'XTickLabel',[]);
            subp2 = subplot(2*(M-1),1,2*m);
            plot(Freqs/1000,imag(squeeze(DC_gen(1,m+1,:))),'-k','LineWidth',2)
            hold on;
            plot(Freqs/1000,imag(squeeze(DC(1,m+1,:))),'-.b','LineWidth',2)
            hold off;
            axis([0 Fs/2000 -1 1]);
            ylabel('Imag')
            if m~=(M-1)
                set(gca,'XTickLabel',[]);
            end
            subp2.Position(2) = subp2.Position(2) + 0.03;
            grid on;
        end
        han=axes(fig,'visible','off');
        han.YLabel.Visible='on';
        han.YLabel.Position(1) = -0.1;
        han.XLabel.Position(2) = -0.05;
        han.XLabel.Visible='on';
        ylabel(han,'Spatial Coherence');
        xlabel(han,'frequency [kHz]');
end