function x = istft(X, hops)
% istft.m computes the inverse Fourier Transform of the STFT matrix X.
%
% Input
%       X:    STFT matrix [Channels x Frequency x Frames]
%       hops: hop size
%
% Output
%       x:    time-domain signal [Samples x Channels]
%
% Author
%       Daniele Mirabilii
%       International Audio Laboratories of Erlangen, Germany
%       email address: daniele.mirabilii@audiolabs-erlangen.de
%
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany

nfft = 2*(length(X(1,:,1))-1); % FFT length (assuming DC and Nyquist are included)

if nargin < 2
    hops = nfft/4; % Default hop size
end

win = hann(nfft); % Hann window (change if different window is desired)

M = size(X,1); % Number of channels
L = size(X,3); % Number of frames
x = zeros(ceil((L-1)*hops+nfft),M); % Pre-allocation time signal

% Frame-wise iFFT
for l = 1:L
    % Symmetric-conjugate spectrum (assuming real-valued signals)
    Xsc = [X(:,:,l).'; conj(X(:,end-1:-1:2,l)).'];

    % Column-wise inverse FFT
    xsc = ifft(Xsc,[],1);

    % Windowing
    xw = win.*xsc(1:nfft,:);

    % OLA index
    idx = 1+(l-1)*hops:nfft+(l-1)*hops;

    % Overlap and add
    x(idx,:) = xw + x(idx,:);
end

% Delete zero-padding which was introduced in stft.m
x = x(1+nfft:end,:);

% Compensate for window gain
gain = sum(win.^2);

% Re-scale signal
x = x.*hops./gain;