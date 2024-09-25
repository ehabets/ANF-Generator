function x = mix_signals(n,C)
% Mix M mutually indepedent signals such that the mixed signals
% exhibit a specific spatial coherence.
%
% Input
%       n : M signals in the STFT domain [Samples x Channels]
%       C : Mixing matrix [Channels x Channels x Frequencies]
%
% Output
%       x : M generated signals [Samples x Channels]
%
% Dependencies
%       stft.m
%       istft.m
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
%------------------------------------------------------------------------------
%
% This code is based on the existing code: https://github.com/ehabets/ANF-Generator
%
% E.A.P. Habets, I. Cohen and S. Gannot, 'Generating nonstationary
% multisensor signals under a spatial coherence constraint,'
% Journal of the Acoustical Society of America, Vol. 124, Issue 5,
% pp. 2911-2917, Nov. 2008.

narginchk(2,2);
M = size(n,2); % Number of sensors
L = size(n,1); % Length input signal
K = (size(C,3)-1)*2; % FFT length
N = permute(stft(n, K),[1 3 2]); % Multi-channel STFT of the input

% Initialization
X = zeros(size(N));  % STFT output matrix

%% Generate output in the STFT domain for each frequency bin k

% Apply mixing matrix
for k = K/2+1:-1:2
    % Filter and sum
    X(:,:,k) = C(:,:,k)' * N(:,:,k);
end

% Inverse STFT
Xp = permute(X,[1 3 2]);
x = real(istft(Xp));
x = x(1:L,:); % Output signals