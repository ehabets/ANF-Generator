function DC = generate_target_coherence(type, params)
% Generate the analytical target coherence based on the model
% (e.g., diffuse spherical), the sensor positions and the
% FFT length. Valid for an arbitrary 3D-array geometry.
%
% Input
%       type:             target coherence model
%                         (spherical, cylindrical, corcos)
%       params.mm:        sensor position matrix [Channels x 3]
%       params.K:         FFT length
%       params.Fs:        sampling frequency
%       params.speed:     wind speed
%       params.direction: wind direction
%
% Output
%       DC:          target coherence matrix [Channels x Channels x K/2+1]
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
% Copyright (c) 2020 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany
% Copyright (c) 2020 Sebastian J. Schlecht

Fs = params.Fs;
K = params.K;
mm = params.mm;
c = 340; % Sound velocity (m/s)

% Angular frequency vector
ww = 2*pi*Fs*(0:K/2)/K;
www = permute(ww, [1 3 2]); % for 3D multiplication

% Matrix of position vectors
rr = permute(mm, [1 3 2]) - permute(mm, [3 1 2]);

% Matrix of inter-sensor distances
d = vecnorm(rr,2,3);

switch type
    case 'corcos'
        % Desired wind speed and direction
        Ud = params.speed; % Wind speed [km/h]
        theta_d = params.direction; % Angle of arrival of the wind stream w.r.t. to the y-axis (North - counterclockwise) in °

        % Speed and direction translations for the Corcos model
        Ud = Ud./3.6;  % Wind speed [m/s]
        U = 0.8*Ud; % Convective turbulence speed
        theta = deg2rad(theta_d); % Angle of arrival in radians
        R = [cos(theta) -sin(theta) 0; ... % Rotation matrix for wind direction
            sin(theta) cos(theta) 0; ...
            0 0 1];
        yy = [ 0 , 1, 0]; % y axis = true North
        u = yy * R;  % Wind direction unit vector
        u_p = [u(2), -u(1), 0]; % Wind direction unit perpendicular vector

        % Coherence parameters
        alpha1 = -0.125;  % Experimental longitudinal decay
        alpha2 = -0.7; % Experimental lateral decay
        dot3 = @(u,rr) sum(permute(u,[1 3 2]).*rr,3); % Vectorized dot product
        alpha__ = alpha1.*(abs(dot3(u,rr)))... % Coherence decay rate all pairs
            + alpha2.*(abs(dot3(u_p,rr)));
        im__ = dot3(u,rr); % Phase difference term all pairs
        AA = (alpha__-1i.*im__).*www;
        DC = exp(AA/U);

    case 'spherical'
        DC = sinc(d.*www/(c*pi));

    case 'cylindrical'
        DC = besselj(0,d.*www/c);
end