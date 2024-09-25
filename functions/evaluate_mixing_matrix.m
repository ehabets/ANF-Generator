function eval = evaluate_mixing_matrix(C_pre,C,decomposition,processing,params,sc_type)
% Analyze the original mixing matrix C_pre (Choleski or Eigenvalue decomp.)
% and the mixing matrix C obtained with one of the proposed processing
% methods.
%
% The function computes:
% 1) Mix balance
% 2) Spectral smoothness
% 3) Coherence error
% 4) IRs (IDFT mixing matrix)
%
% Finally, it plots the performance measures.
%
% Input
%       C_pre         : original mixing matrix [Channels x Channels x Frequencies]
%       C             : processed mixing matrix [Channels x Channels x Frequencies]
%       decomposition : 'CHD' or 'EVD' for Choleski or Eigenvalue Decomposition
%       processing    : 'standard', 'balanced', 'smooth', 'balanced+smooth'
%       params        : spatial coherence parameters
%       sc_type       : type of spatial coherence
%
% Output
%       eval          : struct with all performance measures
%
% Author
%       Daniele Mirabilii
%       International Audio Laboratories of Erlangen, Germany
%       daniele.mirabilii@audiolabs-erlangen.de
%
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany

M = size(C,1); % Number of sensors
K = size(C,3); % FFT length

% Compute frequency-wise balance
bal_pre=zeros(1,K);
bal=zeros(1,K);
for k = 1:K
    bal_pre(k) = sum(abs(C_pre(:,:,k)),'all')/(M*sqrt(M));
    bal(k) = sum(abs(C(:,:,k)),'all')/(M*sqrt(M));
end

% Compute mean balance in dB
bal_pre_mean = mag2db(mean(bal_pre));
bal_mean = mag2db(mean(bal));

% Compute frequency-wise smoothness
smooth_pre = sum(sum(abs(diff(C_pre, 1, 3 )).^2, 1), 2);
smooth = sum(sum(abs(diff(C, 1, 3 )).^2, 1), 2);

% Compute mean smoothness in dB
smooth_pre_mean = pow2db(mean(smooth_pre));
smooth_mean = pow2db(mean(smooth));

% Compute coherence error
stringdm = strcat([decomposition,' ', processing]);
[xi_pre, xi, c_pre, c] = coherence_error(C_pre,C,params,sc_type,stringdm);

% Output
eval.balance_pre = bal_pre_mean; % Balance of C_pre
eval.balance = bal_mean; % Balance of C
eval.smoothness_pre = smooth_pre_mean; % Smoothness of C_pre
eval.smoothness = smooth_mean; % Smoothness of C
eval.coherror_pre = xi_pre; % Coherence error of C_pre
eval.coherror = xi; % Coherence error of C
eval.ir_pre = c_pre; % IRs of C_pre
eval.ir = c; % IR of C

% Summary of output
fprintf('\nImprovements:\n');
fprintf('Spectral Variation: %2.1f dB\n',smooth_mean - smooth_pre_mean);
fprintf('Mix Balance: %2.1f dB\n',bal_mean - bal_pre_mean);
fprintf('Coherence Error: %2.1f dB\n',xi - xi_pre);

% Choose which IRs of the mixing matrix to plot
row = 2;
column = 2;

% Normalize IRs
c_pre = squeeze(abs(c_pre(row,column,:))/max(abs(c_pre(row,column,:))));
c = squeeze(abs(c(row,column,:))/max(abs(c(row,column,:))));

% Plot Spectral Smoothness and Mix Balance
figure()
subplot(2,1,1)
plot(squeeze(smooth),'LineWidth',2); hold on
plot(squeeze(smooth_pre),'-.','LineWidth',2); hold off
legend(sprintf('%s %s: \n \\epsilon = %2.1f dB, \\beta = %2.1f dB', decomposition, processing, smooth_mean, bal_mean),...
    sprintf('%s: \n \\epsilon = %2.1f dB, \\beta = %2.1f dB',decomposition, smooth_pre_mean, bal_pre_mean))
grid on;
axis([1 K/2+1 0 max(max([smooth, smooth_pre]))]);
ylabel('Smoothness');
set(gca,'XTickLabel',[]);
subplot(2,1,2)
plot(bal, 'LineWidth',2); hold on,
plot(bal_pre,'-.', 'LineWidth',2); hold off,
grid on;
xlabel('frequency');
ylabel('Balance');
axis([1 K/2+1 0 1]);

% Plot IRs of the mixing matrix
figure()
subplot(1,2,1)
plot(mag2db(c_pre), 'LineWidth', 1);
title(sprintf('%s - IR(%d,%d)',decomposition, row, column)), grid on, axis tight
xlabel('Time');
ylabel('Magnitude [dB]');
axis([1 2*(K-1)  min(min([mag2db(c_pre), mag2db(c)])) 0]);
subplot(1,2,2)
plot(mag2db(c), 'LineWidth', 1);
title(sprintf('%s %s - IR(%d,%d)',decomposition, processing, row, column)), grid on, axis tight
xlabel('Time');
ylabel('Magnitude [dB]');
axis([1 2*(K-1)  min(min([mag2db(c_pre), mag2db(c)])) 0]);

function [xi_pre, xi, c_pre, c] = coherence_error(C_pre,C,params,sc_type,stringdm)
% Compute the coherence error yielded by the mixing matrix C using the higher
% DFT length K1, as the MSE between the target coherence DC1 and C1'*C1.
%
% Dependencies: generate_target_coherence.m
%
% Input
%        C:        post-processed mixing matrix (smoothed, balanced or both)
%        C_pre:    pre-processed mixing matrix
%        params:    spatial coherence parameters
%        sc_type:  type of spatial coherence
%        stringdm: string for the plots (name of the mixing matrix, e.g., CHD)
% Output
%        nMSE_pre: coherence error pre-processed mixing matrix
%        nMSE:     coherence error post-processed mixing matrix
%        c_pre:    IR matrix of the pre-processed mixing matrix
%        c:        IR matrix of the post-processed mixing matrix

K = 2*(length(C(1,1,:))-1); % Current DFT length
M = length(C(:,1,1)); % Number of channels
K1 = floor(pi*K); % Increased DFT length K1>K

% Pre-allocate memory
cc_pre = zeros(M,M,K);
cc = zeros(M,M,K);
c_pre = zeros(M,M,K);
c = zeros(M,M,K);
C1_pre = zeros(M,M,K1);
C1 = zeros(M,M,K1);

for p = 1:M
    for q = 1:M
        % Conjugate symmetric of the spectrum
        C_pre(p,q,K/2+2:K) = conj(C_pre(p,q,K/2:-1:2));
        C(p,q,K/2+2:K) = conj(C(p,q,K/2:-1:2));

        % IFFT
        cc_pre(p,q,:) = ifft(C_pre(p,q,:));
        cc(p,q,:) = ifft(C(p,q,:));

        % Circular shift
        c_pre(p,q,:) = circshift(squeeze(cc_pre(p,q,:)),K/2);
        c(p,q,:) = circshift(squeeze(cc(p,q,:)),K/2);

        % FFT with increased frame length K1
        C1_pre(p,q,:) = fft(squeeze(c_pre(p,q,:)), K1);
        C1(p,q,:) = fft(squeeze(c(p,q,:)), K1);
    end
end

% Generate target spatial coherence with DFT length K1>K
params.K = K1;
DC1 = generate_target_coherence(sc_type, params);

% Compute generated coherence with DFT-length K1 as C'*C and nMSE
G_pre = zeros(M,M,K);
G = zeros(M,M,K);
nMSEk_pre = zeros(1,K);
nMSEk = zeros(1,K);

for k = 1:K1/2+1
    G_pre(:,:,k) = C1_pre(:,:,k)'*C1_pre(:,:,k);
    G(:,:,k) = C1(:,:,k)'*C1(:,:,k);
    nMSEk_pre(k) = norm(G_pre(:,:,k) - DC1(:,:,k), 'fro').^2;
    nMSEk(k) = norm(G(:,:,k) - DC1(:,:,k), 'fro').^2;
end

% Compute normalized mean squared error between target and generated coherence
xi_pre = pow2db(mean(nMSEk_pre));
xi = pow2db(mean(nMSEk));

% Plot example of pair-wise coherence iDFT - DFT with K1>K (choose desired pair)
Freqs = linspace(0,params.Fs/2,K1/2+1); %vector of frequencies in Hz
decomposition = stringdm(1:3);
processing = stringdm(5:end);
if strcmp(sc_type, 'spherical') || strcmp(sc_type, 'cylindrical')
    % Case diffuse (real-valued coherence)
    fig = figure();
    subplot(2,1,1);
    plot(Freqs/1000,real(squeeze(G_pre(1,2,:))),'-k','LineWidth',2)
    hold on;
    plot(Freqs/1000,real(squeeze(DC1(1,2,:))),'-.b','LineWidth',2)
    hold off;
    sgtitle('Coherence error - Sensors 1-2');
    title(sprintf('%s',decomposition));
    axis([0 params.Fs/2000 -1 1]);
    grid on;
    legend(sprintf('\\xi = %2.1f dB', xi_pre),'Target');
    set(gca,'XTickLabel',[]);
    subplot(2,1,2);
    plot(Freqs/1000,real(squeeze(G(1,2,:))),'-k','LineWidth',2)
    hold on;
    plot(Freqs/1000,real(squeeze(DC1(1,2,:))),'-.b','LineWidth',2)
    hold off;
    title(sprintf('%s %s',decomposition, processing));
    axis([0 params.Fs/2000 -1 1]);
    legend(sprintf('\\xi = %2.1f dB', xi),'Target');
    grid on;
    han=axes(fig,'visible','off');
    han.YLabel.Visible='on';
    han.YLabel.Position(1) = -0.1;
    han.XLabel.Position(2) = -0.05;
    han.XLabel.Visible='on';
    ylabel(han,'Spatial coherence');
    xlabel(han,'Frequency [kHz]');
else
    % Case 'corcos' (or generally complex-valued coherence)
    fig = figure();
    subplot(2,2,1)
    plot(Freqs/1000,real(squeeze(G_pre(1,2,:))),'-k','LineWidth',2);
    hold on;
    plot(Freqs/1000,real(squeeze(DC1(1,2,:))),'-.b','LineWidth',2)
    hold off;
    grid on;
    axis([0 params.Fs/2000  -1 1]);
    ylabel('Real')
    set(gca,'XTickLabel',[]);
    title(sprintf('%s',decomposition));
    legend(sprintf('\\xi = %2.1f dB', xi_pre),'Target');
    set(gca,'XTickLabel',[]);
    subplot(2,2,3)
    plot(Freqs/1000,imag(squeeze(G_pre(1,2,:))),'-k','LineWidth',2)
    hold on;
    plot(Freqs/1000,imag(squeeze(DC1(1,2,:))),'-.b','LineWidth',2)
    hold off;
    axis([0 params.Fs/2000 -1 1]);
    ylabel('Imag')
    grid on;
    subplot(2,2,2)
    plot(Freqs/1000,real(squeeze(G(1,2,:))),'-k','LineWidth',2);
    hold on;
    plot(Freqs/1000,real(squeeze(DC1(1,2,:))),'-.b','LineWidth',2)
    hold off;
    grid on;
    axis([0 params.Fs/2000  -1 1]);
    set(gca,'XTickLabel',[]);
    title(sprintf('%s %s',decomposition,processing));
    legend(sprintf('\\xi = %2.1f dB', xi),'Target');
    set(gca,'XTickLabel',[]);
    subplot(2,2,4)
    plot(Freqs/1000,imag(squeeze(G(1,2,:))),'-k','LineWidth',2)
    hold on;
    plot(Freqs/1000,imag(squeeze(DC1(1,2,:))),'-.b','LineWidth',2)
    hold off;
    axis([0 params.Fs/2000 -1 1]);
    grid on;
    han=axes(fig,'visible','off');
    sgtitle('Coherence error - Sensors 1-2');
    han.YLabel.Visible='on';
    han.YLabel.Position(1) = -0.1;
    han.XLabel.Position(2) = -0.05;
    han.XLabel.Visible='on';
    ylabel(han,'Spatial Coherence');
    xlabel(han,'Frequency [kHz]');
end