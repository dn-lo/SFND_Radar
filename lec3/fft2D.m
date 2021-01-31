% 2-D Transform
% The 2-D Fourier transform is useful for processing 2-D signals and other 2-D data such as images.
% Create and plot 2-D data with repeated blocks.

P = peaks(20);
X = repmat(P,[5 10]);
imagesc(X)

% Compute the 2-D Fourier transform of the data.  
% Shift the zero-frequency component to the center of the output, and 
% plot the resulting 100-by-200 matrix, which is the same size as X.

% Run the 2D FFT across both the dimensions.
P_fft = fft2(P);

% Shift zero-frequency terms to the center of the array
P_fft = fftshift(P_fft);

% Take the absolute value
P_fft = abs(P_fft);

% it can be plotted as an image. Hence, we use the imagesc 
imagesc(P_fft);


