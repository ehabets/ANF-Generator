function B = balanced(C, type)
% Compute optimal balanced B = U * C with max l1( U * C ) for U unitary
%
% Input
%         C: mixing matrix [Channels x Channels]
%         type: 'unitary' (default) for U unitary, and 'orthogonal' for U orthogonal
% Output
%         B: balanced mixing matrix [Channels x Channels]
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

if nargin == 1
    type = 'unitary';
end

N = size(C,1); % Number of channels

% Init with complex Hadamard matrix
% H = dftmtx(N)/sqrt(N);

% Init with a matrix of scaled ones
H = ones(N)/sqrt(N);

% Perform optimization
switch type
    case 'unitary'
        U = absoluteUnitaryProcrustes(C, H);
    case 'orthogonal'
        U = absoluteOrthogonalProcrustes(C, H);
end

% Compute balanced matrix
B = U * C;

function [bestMatrix,l1,U] = absoluteOrthogonalProcrustes(A, B)
% Solves || abs(U * A) - abs(B) ||_F with orthogonal U for A, B and U being
% real-valued

MaximumTrails = 20; % Number of trials
N = size(A,1); % Number of channels
newSigns = sign(real(B)+eps); % Sign matrix of B

U = cell(1,MaximumTrails);
l1 = zeros(1,MaximumTrails);

for it = 1:MaximumTrails
    % Perform optimization
    U{it} = real(signVariableExchange(A, newSigns));

    % Evaluate l1-norm of obtained unitary matrix
    l1(it) = l1norm(U{it}*A);

    % Generate new (random) normalized phase matrix
    newSigns = sign(rand(N)*2-1);
end

% Search matrix (U*A) with max l1-norm
[~,ind] = max(l1);
bestMatrix = U{ind};


function [bestMatrix,l1,U] = absoluteUnitaryProcrustes(A, B)
% Solves || abs(U * A) - abs(B) ||_F with unitary U

MaximumTrails = 20; % Number of trials
N = size(A,1); % Number of channels
newPhases = exp(1i*angle(B)); % Phase matrix of B

U = cell(1,MaximumTrails);
l1 = zeros(1,MaximumTrails);

for it = 1:MaximumTrails
    % Perform optimization
    U{it} = signVariableExchange(A, newPhases);

    % Evaluate l1-norm of obtained unitary matrix
    l1(it) = l1norm(U{it}*A);

    % Generate new (random) normalized phase matrix
    newPhases = exp(1i*rand(N)*2*pi);
end

% Search matrix (U*A) with max l1-norm
[~,ind] = max(l1);
bestMatrix = U{ind};


function [U, l1Fit, Us] = signVariableExchange(A, B)
% Perform the variable exchange optimization

maxIter = 1000; % Max number of iterations
bestFit = Inf; % Best fit in the Frobenius-norm sense
Phases = exp(1i*angle(B)); % Phase matrix of B

newFit = zeros(1,maxIter);
l1Fit = zeros(1,maxIter);
Us = cell(1,maxIter);
for counter = 1:maxIter
    % New phase matrix
    phaseB = Phases;

    % Procrustes solution
    U = procrustes(A,phaseB);

    % Compute Frobenius norm of (U*A-phaseB)
    newFit(counter) = norm(U*A - phaseB,'fro');

    % Check convergence
    improvement = newFit(counter) - bestFit;

    % Update bestFit
    bestFit = newFit(counter);

    % Store l1norm
    l1Fit(counter) = l1norm(U*A);

    % Store the unitary matrix
    Us{counter} = U;

    % Convergence
    if improvement > sqrt(eps)
        warning('Not Improved');
    elseif improvement > -10^-5 % Optional message after converging
        fprintf('Variable exchange optimization converged within %d iterations\n', counter)
        break;
    end

    % Update phase matrix
    Phases = exp(1i*angle(U*A));

    if(counter == maxIter) % optional message for not converging
        warning('Variable exchange optimization did not converge within %d iterations.', maxIter)
    end
end


function n = l1norm(X)
% Normalized l1 norm (element-wise)
% For ||X||_F = 1, the norm is in [0,1].

N = size(X,1); % Number of channels

% Normalized l1-norm
n = sum(abs(X),'all') / (N*sqrt(N));