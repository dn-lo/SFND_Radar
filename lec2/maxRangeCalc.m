% Operating frequency (Hz)
fc = 77.0e9;

% Transmitted power (W)
Ps = 3e-3;      % +5  dBm

% Antenna Gain (linear gain)
G =  10000;     % +40 dBi

% Minimum Detectable Power (W)
Pe = 1e-10;     % -70 dBm

% RCS of a car (m^2)
RCS = 100;      % +20 dB

% Speed of light
c = 3e8;

% Calculate the wavelength
lambda = c/fc;

% Measure the Maximum Range a Radar can see. 
R = (Ps * G^2 * lambda^2 * RCS / (Pe * (4*pi)^3))^(1/4)