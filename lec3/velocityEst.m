%% Calculate the range of four targets with measured beat frequencies
fb_targets = [0, 1.1, 13, 24] * 1e6; % (Hz)

%% Parameters
R_max = 300;     % (m): max radar range
d_res = 1;      % (m): radar range resolution
c = 3e8;        % (m/s): speed of light

%% Calculations

% Find the sweep Bandwidth of the chirp for 1 m resolution
B_sweep = c / (2 *  d_res); % (Hz)

% Calculate the chirp time based on the Radar's Max Range
Ts = 5.5 * 2 * R_max / c; % (s)

% Define the frequency shifts 
calculated_range = c * Ts * fb_targets / (2 * B_sweep); % (m)

% Display the calculated range
disp(calculated_range);