function Cbs = balanced_preserving_smoothness(C)
% Compute balanced Cbs = U * C with high l1( U * C ) for U UNITARY while
% preserving smoothness (=single iteration of balanced.m with closest
% phases of C)
%
% Input
%       C: mixing matrix [Channels x Channels]
%
% Output
%       B: smoothness-preserving balanced  mixing matrix [Channels x Channels]
%
% Dependencies
%       procrustes.m
%
% Authors
%       Sebastian Jiro Schlecht
%       Aalto University, Finland
%       sebastian.schlecht@aalto.fi
%
%       Daniele Mirabilii
%       International Audio Laboratories of Erlangen, Germany
%       daniele.mirabilii@audiolabs-erlangen.de
%
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany
% Copyright (c) 2020 Sebastian J. Schlecht

% Closest phases (phases of C)
B = exp(1i*angle(C));

% Procrustes solution
U = procrustes(C,B);

% Smoothing-preserving balanced matrix
Cbs = U*C;