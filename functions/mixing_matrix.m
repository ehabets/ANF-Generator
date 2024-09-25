function [C, C_none] = mixing_matrix(DC,decomposition,processing)
% Compute the mixing matrix is obtained by decomposing the spatial coherence
% matrix using the Choleski (CHD) or the Eigenvalue (EVD) decomposition.
% Optionally, enhances different properties of the obtained mixing matrix
% by selecting the dedicated processing method, i.e., 'smooth' enhances the
% spectral smoothness, 'balanced' enhnaces the balance and 'balanced+smooth'
% enhances both jointly. The processing method 'standard' leaves the mixing
% matrix as obtained by the CHD or EVD unaltered.
%
% Input
%       DC            : Desired coherence [Channels x Channels x Frequencies]
%       decomposition : 'CHD' or 'EVD' for Choleski or Eigenvalue Decomposition
%       processing    : 'standard', 'balanced', 'smooth', 'balanced+smooth'
%
% Output
%       C             : Mixing matrix [Channels x Channels x Frequencies]
%
% Dependencies
%       smoothness.m
%       balanced.m
%       balanced_preserving_smoothness.m
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

narginchk(1,3);
if nargin < 2
    decomposition = 'CHD';
    processing = 'balanced+smooth';
end
if nargin < 3
    processing = 'balanced+smooth';
end

% Initialization
M = size(DC,1); % Number of sensors
K = (size(DC,3)-1)*2; % FFT length
C = zeros(size(DC)); % STFT mixing matrix

% Direct current component definition for the mixing matrix (fully coherent)
C(:,:,1) = ones(M,M)./sqrt(M);

%% Generate output in the STFT domain for each frequency bin k

% Apply CHD or EVD
for k = K/2+1:-1:2
    switch lower(decomposition)
        case 'chd'
            if rank(DC(:,:,k)) == M
                C(:,:,k) = chol(DC(:,:,k)); % Cholesky decomposition
            else
                % Fallback solution in case the Cholesky decomposition cannot be computed
                [V,D] = eig(DC(:,:,k)); % Eigenvalue decomposition
                C(:,:,k) = sqrt(D) * V';
            end
        case 'evd'
            [V,D] = eig(DC(:,:,k)); % Eigenvalue decomposition
            C(:,:,k) = sqrt(D) * V';
        otherwise
            error('Unknown decomposition method specified. Please select "CHD" or "EVD".');
    end
end

% Store original CHD/EVD mixing matrix for evaluation
C_none = C;

% Apply selected processing method
for k = K/2+1:-1:2
    switch processing
        case 'standard'
            continue % No post-processing
        case 'balanced'
            if k == K/2+1
                C(:,:,k) = balanced(C(:,:,k),'orthogonal'); % Balance initial state
            else
                C(:,:,k) = balanced(C(:,:,k),'unitary'); % Induce balance
            end
        case 'smooth'
            C(:,:,k-1) = smoothness(C(:,:,k),C(:,:,k-1)); % Induce smoothness
        case 'balanced+smooth'
            if k == K/2+1
                C(:,:,k) = balanced(C(:,:,k),'orthogonal'); % Balance initial state
            end
            C(:,:,k-1) = smoothness(C(:,:,k),C(:,:,k-1)); % Induce smoothness
            C(:,:,k-1) = balanced_preserving_smoothness(C(:,:,k-1)); % Induce smoothness-preserving balance
        otherwise
            error('Unknown processing method specified. Please select "standard", "balanced", "smooth", or "balanced+smooth."');
    end
end