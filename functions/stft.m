function X = stft(x, nfft, hops)
% stft.m computes the Short Time Fourier Transform of a multi-channel input
%
% Input
%       x:       time signal [Sample x Channels]
%       nft:     FFT length
%       win_len: window length
%       hops:    hop size
%
% Output
%       X:       STFT matrix [Channels x Frequencies x Frames]
%
% Note
%       If the input signal is mono-channel, the STFT matrix will have size
%       [1 x Frequencies x Frames]. Use squeeze(X) to remove the singleton
%       dimension.
%
% Author
%       Daniele Mirabilii
%       International Audio Laboratories of Erlangen, Germany
%       daniele.mirabilii@audiolabs-erlangen.de
%
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany

if nargin < 2
    nfft = 1024;
end
if nargin < 3
    hops = nfft/4;
end
win = hann(nfft); % Hann window (change if different window is desired)
M = size(x,2); % Number of channels

% Delay to avoid signal smearing due to the windowing at the beginning/end
xd = [zeros(nfft,M); x; zeros(nfft,M)];
T = size(xd,1); % Number of samples including the introduced delay

K = floor(nfft/2)+1; % Number of frequency instances
L = floor((T-nfft)/hops)+1; % Number of frames

% Fit the delayed signal into the actual number of frames
% xc = xd;%xd(1:(L-1)*hops+nfft,:);

% Pre-allocation STFT matrix
X = zeros(M,K,L);

% Frame-wise FFT
for l = 0:L-1
    % Frame-wise index
    idx = 1+l*hops:nfft+l*hops;

    % Windowing
    xw = win.*xd(idx,:);

    % Column-wise FFT of the windowed-signal matrix
    Xw = fft(xw,nfft,1);

    % Allocation in the STFT matrix
    X(:,:,l+1) = Xw(1:K,:).';
end