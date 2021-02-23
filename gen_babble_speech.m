% DESCRIPTION:
%   Example that generates babble speech received by a uniform linear
%   array of sensors.
%
% REFERENCE:
%   E.A.P. Habets, I. Cohen and S. Gannot, 'Generating
%   nonstationary multisensor signals under a spatial
%   coherence constraint', Journal of the Acoustical Society
%   of America, Vol. 124, Issue 5, pp. 2911-2917, Nov. 2008.
%
% REMARKS:
%   For babble speech the Cholesky decomposition is preferred over the
%   eigenvalue decomposition.
%
% REVISION HISTORY:
%   2008 - E.A.P. Habets
%       * Initial implementation
%   25/01/2021 - E.A.P. Habets
%       * Removed third-party dependencies

close all;
clear;

% Initialization
Fs = 8000; % Sample frequency (Hz)
c = 340; % Sound velocity (m/s)
K = 256; % FFT length
M = 4; % Number of sensors
d = 0.1; % Inter sensor distance (m)
type_nf = 'spherical'; % Type of noise field: 'spherical' or 'cylindrical'
L = 20*Fs; % Data length

%% Generate M mutually 'independent' babble speech input signals
[data,Fs_data] = audioread('babble_8kHz.wav');
if Fs ~= Fs_data
    error('Sample frequency of input file is incorrect.');
end
data = data - mean(data);
babble = zeros(L,M);
for m=1:M
    babble(:,m) = data((m-1)*L+1:m*L);
end

%% Generate matrix with desired spatial coherence
ww = 2*pi*Fs*(0:K/2)/K;
DC = zeros(M,M,K/2+1);
for p = 1:M
    for q = 1:M
        if p == q
            DC(p,q,:) = ones(1,1,K/2+1);
        else
            switch lower(type_nf)
                case 'spherical'
                    DC(p,q,:) = sinc(ww*abs(p-q)*d/(c*pi));
                    
                case 'cylindrical'
                    DC(p,q,:) = bessel(0,ww*abs(p-q)*d/c);
                    
                otherwise
                    error('Unknown noise field.')
            end
        end
    end
end

if M>=6
    EPS = 1e-10;
    for k = 1:K/2+1
        DC(:,:,k) = DC(:,:,k) + EPS*eye(M);
    end
end

%% Generate sensor signals with desired spatial coherence
x = mix_signals(babble,DC,'cholesky');

%% Compare desired and generated coherence
K_eval = 256;
ww = 2*pi*Fs*(0:K_eval/2)/K_eval;
sc_theory = zeros(M-1,K_eval/2+1);
sc_generated = zeros(M-1,K_eval/2+1);

% Calculalte STFT and PSD of all output signals
X = stft(x,'Window',hanning(K_eval),'OverlapLength',0.75*K_eval,'FFTLength',K_eval,'Centered',false);
X = X(1:K_eval/2+1,:,:);
phi_x = mean(abs(X).^2,2);

% Calculate spatial coherence of desired and generated signals
for m = 1:M-1
    switch lower(type_nf)
        case 'spherical'
            sc_theory(m,:) = sinc(ww*m*d/(c*pi));
            
        case 'cylindrical'
            sc_theory(m,:) = bessel(0,ww*m*d/c);
    end
    
    % Compute cross-PSD of x_1 and x_(m+1)
    psi_x =  mean(X(:,:,1) .* conj(X(:,:,m+1)),2);
    
    % Compute real-part of complex coherence between x_1 and x_(m+1)
    sc_generated(m,:) = real(psi_x ./ sqrt(phi_x(:,1,1) .* phi_x(:,1,m+1))).';
end

% Calculate normalized mean square error
NMSE = zeros(M,1);
for m = 1:M-1
    NMSE(m) = 10*log10(sum(((sc_theory(m,:))-(sc_generated(m,:))).^2)./sum((sc_theory(m,:)).^2));
end

% Plot spatial coherence of two sensor pairs
figure(1);
MM=min(2,M-1);
Freqs=0:(Fs/2)/(K/2):Fs/2;
for m = 1:MM
    subplot(MM,1,m);
    plot(Freqs/1000,sc_theory(m,:),'-k','LineWidth',1.5)
    hold on;
    plot(Freqs/1000,sc_generated(m,:),'-.b','LineWidth',1.5)
    hold off;
    xlabel('Frequency [kHz]');
    ylabel('Real(Spatial Coherence)');
    title(sprintf('Inter sensor distance %1.2f m',m*d));
    legend('Theory',sprintf('Proposed Method (NMSE = %2.1f dB)',NMSE(m)));
    grid on;
end

% Save babble speech
audiowrite('mc_babble_speech_example.wav',x,Fs);
