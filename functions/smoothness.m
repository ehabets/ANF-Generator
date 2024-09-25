function C_n = smoothness(C, C_plus)
% Compute nearest UNITARY matrix U (in the Frobenius norm) via Procrustes
% solution such that (U*C_plus - C) is minimized. The output is the updated
% matrix C_n = U*C_plus.
%
% Input
%       C: mixing matrix [Channels x Channels]
%
% Output
%       B: smooth (w.r.t. neighbour C) mixing matrix [Channels x Channels]
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

U = procrustes(C_plus,C);
C_n = U*C_plus;