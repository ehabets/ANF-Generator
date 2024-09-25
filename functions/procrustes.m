function U = procrustes(A,B)
% Compute nearest orthogonal matrix U (in the Frobenius norm) via SVD such
% that it minimizes U*A - B
%
% Input
%        A: first matrix
%        B: second matrix
%
% Output
%        U: unitary matrix that minimizes (U*A - B) in the Frobenius norm
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

% Compute SVD of BA'
[W,~,V] = svd(B*A');

% Compute orthogonal polar factor of (BA')
U = W*V';