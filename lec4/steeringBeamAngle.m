% Operating frequency (Hz)
fc = 77.0e9;

% Phase increment (deg)
phi = 45;

% Speed of light (m/s)
c = 3e8;

% Calculate the wavelength (m)
lambda = c/fc;

% Antenna element spacing (m)
d = lambda/2;

% Steering angle of antenna beam (deg)
theta = asind(phi/360 * lambda/d)
