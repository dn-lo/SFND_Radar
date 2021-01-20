%% Use a Fourier transform to find the frequency components of a signal buried in noise
%% Parameters
Fs = 1000;            % (Hz): Sampling frequency                    
T = 1/Fs;             % (s): Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

%% Define signal characteristics
% Signal
S = 0.7 * sin(2*pi * 77 * t) + 2 * sin(2*pi * 43 * t);
% Corrupt with noise
X = S + 2 * randn(size(t));

% Plot the noisy signal in the time domain. It is difficult to identify the frequency components by looking at the signal X(t). 
plot(1000*t(1:50) ,X(1:50), 'linewidth', 1.5)
grid on
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (ms)')
ylabel('X(t)')

%% Compute the spectrum of the signal with FFT

% Extract the L-point DFT of the signal (complex)
signal_fft = fft(X, L);

% Compute the double-sided spectrum P2 (real): corresponds to the 
% (scaled) magnitude of the FFT
P2 = abs(signal_fft/L);

% Compute single-sided spectrum of the signal
P1 = P2(1:L/2+1);
% P1(:, 2:end-1) = 2 * P1(:, 2:end-1);

% Plotting
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
