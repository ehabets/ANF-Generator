function x = mix_signals(n,DC,method)
% SYNTAX:
%   x = mix_signals(n,DC,method)
%
% DESCRIPTION:
%   Mix M mutually indepedent signals n of length L such that the
%   mixed signals exhibit a specific spatial coherence given by DC.
%
% INPUT:
%   x - (double) M signals in the time domain [L x M]
%   DC - (double) Desired spatial coherence [M x M x K/2+1]
%   method - 'cholesky' (default) or 'eigen'
%
% OUTPUT:
%   x - (double) M output signals of length L in the time domain [L x M]
%
% ASSUMPTIONS AND LIMITATIONS:
%   Assumes output signals are fully coherent at 0 Hz
%
% REFERENCE:
%   E.A.P. Habets, I. Cohen and S. Gannot, 'Generating
%   nonstationary multisensor signals under a spatial
%   coherence constraint', Journal of the Acoustical Society
%   of America, Vol. 124, Issue 5, pp. 2911-2917, Nov. 2008.
%
% REVISION HISTORY:
%   2008 - E.A.P. Habets
%       * Initial implementation
%   25/01/2021 - E.A.P. Habets
%       * Changed license to MIT
%       * Removed third-party dependencies

% MIT License
%
% Copyright (C) 2021 E.A.P. Habets
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

narginchk(2,3);

if nargin < 3
    method = 'cholesky';
end

M = size(n,2); % Number of sensors
K = (size(DC,3)-1)*2; % Number of frequency bins

% Compute short-time Fourier transform (STFT) of all input signals
n = [zeros(K/2,M) ; n ; zeros(K/2,M)];
N = stft(n,'Window',hanning(K),'OverlapLength',0.75*K,'FFTLength',K,'Centered',false);

% Generate output signal in the STFT domain for each frequency bin k
X = zeros(size(N)); % STFT output matrix
for k = 2:K/2+1
    switch lower(method)
        case 'cholesky'
            C = chol(DC(:,:,k));
            
        case 'eigen'
            [V,D] = eig(DC(:,:,k));
            C = sqrt(D) * V';
            
        otherwise
            error('Unknown method specified.');
    end
    
    X(k,:,:) = squeeze(N(k,:,:)) * conj(C);
end
X(K/2+2:end,:,:) = conj(X(end:-1:K/2+2,:,:));

% Compute inverse STFT
x = istft(X,'Window',hanning(K),'OverlapLength',0.75*K,'FFTLength',K,'Centered',false,'ConjugateSymmetric',true);
x = x(K/2+1:end-K/2,:);