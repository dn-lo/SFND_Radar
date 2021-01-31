%% Calculate the range of four targets with measured doppler frequency shifts
fd_targets = [3, -4.5, 11, -3] * 1e3; % (Hz)

%% Parameters
f = 77e9;       % (Hz): radar operating frequency
c = 3e8;        % (m/s): speed of light

%% Calculations

% Calculate the wavelength
lambda = c/f;   % (m)

% Calculate the velocity of the targets  fd = 2*vr/lambda
vr = fd_targets * lambda / 2;

% % Display results
% disp(lambda)
% disp(vr)

%% Question
dist0 = 200;    % (m): initial targets distance
v_ego = 5;      % (m/s): ego vehicle speed
t = 5;          % (s): time

% Current target distances
dist = dist0 - (v_ego + vr) * t; %(m)
disp(dist)

