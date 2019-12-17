function [Cxy, f] = complex_cohere(varargin)
%COMPLEX_COHERE Complex coherence function estimate.
%   Cxy = COMPLEX_COHERE(X,Y,NFFT,Fs,WINDOW,NOVERLAP) estimates the 
%   complex coherence of X and Y using Welch's averaged periodogram method.  
%   The complex coherence is a function of frequency with values between 
%   0 and 1 that indicate how well the input X corresponds to the output Y
%   at each frequency.  X and Y are  divided into overlapping sections by
%   NOVERLAP samples, each of which is detrended, then windowed by the 
%   WINDOW parameter, then zero-padded to length NFFT.  
%   The magnitude squared of the length NFFT DFTs of the sections of X and 
%   the sections of Y are averaged to form Pxx and Pyy, the Power Spectral
%   Densities of X and Y respectively. The products of the length NFFT DFTs
%   of the sections of X and Y are averaged to form Pxy, the Cross Spectral
%   Density of X and Y. The coherence Cxy is given by
%       Cxy = Pxy./sqrt(Pxx.*Pyy)
%   Cxy has length NFFT/2+1 for NFFT even, (NFFT+1)/2 for NFFT odd, or NFFT
%   if X or Y is complex. Fs is the sampling frequency which does
%   not effect the cross spectrum estimate but is used for scaling of plots.
%
%   COMPLEX_COHERE with no output arguments plots the coherence in the current 
%   figure window.

narginchk(6,6);
x = varargin{1};
y = varargin{2};
nfft = varargin{3};
Fs = varargin{4};
window=varargin{5};
noverlap = varargin{6};

% compute PSD and CSD
window = window(:);
n = length(x);		% Number of data points
nwind = length(window); % length of window
if n < nwind    % zero-pad x , y if length is less than the window length
    x(nwind)=0;
    y(nwind)=0;  
    n=nwind;
end
x = x(:);		% Make sure x is a column vector
y = y(:);		% Make sure y is a column vector
k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;

Pxx = zeros(nfft,1); 
Pyy = zeros(nfft,1); 
Pxy = zeros(nfft,1); 
for i=1:k
    xw = window.*x(index);
    yw = window.*y(index);
    index = index + (nwind - noverlap);
    Xx = fft(xw,nfft);
    Yy = fft(yw,nfft);
    Xx2 = abs(Xx).^2;
    Yy2 = abs(Yy).^2;
    Xy2 = Yy.*conj(Xx);
    Pxx = Pxx + Xx2;
    Pyy = Pyy + Yy2;
    Pxy = Pxy + Xy2;
end

% Select first half
if ~any(any(imag([x y])~=0))   % if x and y are not complex
    if rem(nfft,2)    % nfft odd
        select = 1:(nfft+1)/2;
    else
        select = 1:nfft/2+1;   % include DC AND Nyquist
    end
    Pxx = Pxx(select);
    Pyy = Pyy(select);
    Pxy = Pxy(select);
else
    select = 1:nfft;
end
Coh = Pxy./sqrt(Pxx.*Pyy);
freq_vector = (select - 1)'*Fs/nfft;

% set up output parameters
if (nargout == 2)
   Cxy = Coh;
   f = freq_vector;
elseif (nargout == 1)
   Cxy = Coh;
elseif (nargout == 0)   % do a plot
   newplot;
   plot(freq_vector,Coh), grid on
   xlabel('Frequency'), ylabel('Coherence Function Estimate');
end