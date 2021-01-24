function x = mix_signals(n,DC,method)

% Mix M mutually indepedent signals such that the mixed signals 
% exhibit a specific spatial coherence.
%
% Input parameters:
%       n      : M signals in the STFT domain [L x M]
%       DC     : Desired coherence [M x M x K/2+1]
%       method : 'cholesky' or 'eigen'
%
% Output parameters:
%       x      : M generated signals [L x M]
%
% Author       : E.A.P. Habets
% Date         : 29-06-2017
%
% Reference    : E.A.P. Habets, I. Cohen and S. Gannot, 'Generating 
%                nonstationary multisensor signals under a spatial 
%                coherence constraint', Journal of the Acoustical Society
%                of America, Vol. 124, Issue 5, pp. 2911-2917, Nov. 2008.

% MIT License
%
% Copyright (C) 2009-2017 E.A.P. Habets
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
L = size(n,1); % Length input signal
K = (size(DC,3)-1)*2;

% Short-time Fourier transform
for m = 1 : M
    N(m,:,:) = stft(n(:,m), K, K/4, 1).';
end

% Initialization
C = zeros(size(DC)); % STFT mixing matrix
X = zeros(size(N));  % STFT output matrix
X(:,:,1) = X(1,1,1);

% Generate output in the STFT domain for each frequency bin k
for k = 2:K/2+1
    switch lower(method)
        case 'cholesky'
            C(:,:,k) = chol(DC(:,:,k));
            
        case 'eigen'
            [V,D] = eig(DC(:,:,k));
            C(:,:,k) = sqrt(D) * V';
            
        otherwise
            error('Unknown method specified.');
    end
    
    X(:,:,k) = C(:,:,k)' * N(:,:,k);
end

% Inverse STFT
for m = 1 : M
    x(:,m) = real(istft(squeeze(X(m,:,:)).', K, K/4, 1));
end
x = x(1:L,:);