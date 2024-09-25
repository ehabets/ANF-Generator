function CC = mccoherence(x,nfft,hops)
% Compute the complex coherence of a real-signals matrix x using the
% Welch's periodogram approach.
%
% Input
%       x:    multi-channel time signals [Samples x Channels]
%       nfft: FFT length
%       hops: hop size
%
% Output
%       CC:   estimated complex coherence matrix [Channels x Channels x Frequencies]
%
% Dependencies
%       stft.m
%
% Author
%       Daniele Mirabilii
%       International Audio Laboratories of Erlangen, Germany
%       daniele.mirabilii@audiolabs-erlangen.de
%
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany

% Check if channels >1
if length(x(1,:)) == 1
    error('Number of channels must be > 1 to compute the spatial coherence.');
end

if nargin < 2
    nfft = 512; % Default FFT size
end
if nargin < 3
    hops = nfft/4; % Default hops
end

% Compute STFT of the input
X = squeeze(stft(x,nfft,hops));
[M, K, L] = size(X); % Number of channels, frequencies and frames

% Pre-allocation cPSD matrix
Pxx = zeros(K,M,M);

% Frame-wise processing
for l=1:L
    % Compute instantaneous spectrum
    Xf = squeeze(X(:,:,l)).';

    % Compute nstantaneous cross-power spectrum
    XX = Xf .* conj(permute(Xf, [1 3 2]));

    % Update
    Pxx = Pxx + XX;
end

% Complex coherence matrix
CC = zeros(M,M,K);
for r=1:M
    for c=1:M
        CC(r,c,:) = Pxx(:,r,c)./sqrt(Pxx(:,r,r).*Pxx(:,c,c));
    end
end