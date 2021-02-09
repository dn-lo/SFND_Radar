clear all
clc;

%% System requirements
% Max Range (m)
rangeMax = 200;
% Range Resolution (m)
rangeRes = 1;
% Max Velocity (m/s)
velMax = 70;
% Velocity Resolution(m/s)
velRes = 3;
% Speed of light (m/s)
c = 3e8;

%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains contant
% Initial target position (m)
x0 = 155;
% Initial target speed (m/s)
v0 = -30; % receding target
% v0 = 30; % approaching target

%% FMCW Waveform Generation

% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

% Bandwidth (Hz)
B = c / (2 * rangeRes);

% Chirp time (s)
Tchirp = 5.5 * 2 * rangeMax / c;

% Slope (1/s2)
slope = B / Tchirp;

%Operating carrier frequency of Radar (Hz)
fc = 77e9;             %carrier freq
                                      
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd = 128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr = 1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp. 
% We will have a total of Nr*Nd time samples to analyze at each
% update. In this case the last timestamp is t(end) = 1 ms
t = linspace(0, Nd*Tchirp, Nr*Nd+1); %total time for samples
t = t(1:end-1); % cut last point to make it of size Nr*Nd     

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); % transmitted signal
Rx = zeros(1,length(t)); % received signal
Mix = zeros(1,length(t)); % beat signal

%Similar vectors for range_covered and time delay.
r_t = zeros(1,length(t));
td = zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 
for i=1:length(t)         
    
    %For each time stamp update the Range of the Target for constant velocity. 
    %Assuming t0 = t(1) = 0 s
    r_t(i) = x0 + v0 * t(i);
    
    %Compute signal return delay
    td(i) = 2 * r_t(i) / c;

    %For each time sample we need update the transmitted and
    %received signal. NB:  phase of cos is the integral of frequency
    %f=fc+slope*t
    % Assuming sawtooth chirp
    
    % transmitted signal
    Tx(i)  = cos(2*pi * (fc*t(i) + slope*t(i)^2/2));
    
    % received signal
    Rx(i)  = cos(2*pi * (fc*(t(i)-td(i)) + slope*(t(i)-td(i))^2/2));
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal. Check sign of speed
    Mix(i) = Tx(i) .* Rx(i);
    
end

%% RANGE MEASUREMENT


% reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
% Range and Doppler FFT respectively.
Mix_mat = reshape(Mix, [Nr, Nd]);

% run the FFT on the beat signal along the range bins dimension (Nr) and
% normalize.
sig_fft = fft(Mix_mat, Nr, 1);
sig_fft = sig_fft / Nr;

% Take the absolute value of FFT output
sig_fft = abs(sig_fft);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_fft = sig_fft(1:(Nr/2+1));

% Frequency and range vectors
f = 1/(t(2)-t(1)) * (0:(Nr/2))/Nr;
range = c / (2 * slope) * f; 

% plotting the range
figure ('Name','Range from First FFT')

% plot FFT output 
hold on, grid on
plot(range, sig_fft, 'linewidth', 1.2)
 
axis ([0 200 0 0.5]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has response in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Mix values.

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix_mat, Nr, Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2, 1:Nd);
% Shift zero-frequency component
sig_fft2 = fftshift(sig_fft2);
% Obtain Range-Doppler map from FFT modulus and convert it to dB
RDM = abs(sig_fft2);
RDM = pow2db(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100, 100, Nd);
range_axis = linspace(-200 , 200, Nr/2) * ((Nr/2)/400);
figure()
surf(doppler_axis, range_axis, RDM);
colorbar;
xlabel('speed [m/s]')
ylabel('range [m]')

%% CFAR implementation

%Slide Window through the complete Range Doppler Map
%Create empty matrix to hold the result of CA-CFAR 2D
sig_cfar = zeros(size(RDM));

%Select the number of Training Cells in both the dimensions.
Tr = 16;
Td = 8;
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 2;

% offset the threshold by SNR value in dB
SNR = 10;   % dB
threshold_cfar = zeros(size(RDM));

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(size(RDM));

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

% Exclude cells that are too close to the borders: CUT should have margins
% for training and guard cells on both range and Doppler axis
for i = (1+(Tr+Gr)):(size(RDM,1)-(Tr+Gr))
    
    for j = (1+(Td+Gd)):(size(RDM,2)-(Td+Gd))
    % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
    % CFAR
    
    % Here we use linear indexing to select matrix indices
    % Select the grid that includes the training, guard and test cells
    [ii_tot, jj_tot] = meshgrid((i-(Tr+Gr)):(i+(Tr+Gr)), (j-(Td+Gd)):(j+(Td+Gd)));
    ids_tot = sub2ind(size(RDM), ii_tot(:), jj_tot(:));
    % Select the grid that includes cells in the guard region and cell under test.
    [ii_grd, jj_grd] = meshgrid((i-Gr):(i+Gr), (j-Gd):(j+Gd));
    ids_grd = sub2ind(size(RDM), ii_grd(:), jj_grd(:));
    % Subtract sets to obtain the training cells
    ids_trn = setdiff(ids_tot, ids_grd);
    
    % Average noise across training cells (convert to power for averaging)
    % noise is in dB
    noise_level(i,j) = mean(db2pow(RDM(ids_trn)));
    noise_level(i,j) = pow2db(noise_level(i,j));
    
    % Add the offset (in dB) to the threshold
    threshold_cfar(i,j) = noise_level(i,j) + SNR;
    
    % If the CUT signal level is greater than the Threshold, assign a value 
    % of 1, else equate it to zero.
    sig_cfar(i,j) = double( RDM(i,j) > threshold_cfar(i,j) );
    
    end
end

% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
% Already done by initialization


%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure()
surf(doppler_axis, range_axis, sig_cfar)
colorbar;
xlabel('speed [m/s]')
ylabel('range [m]')


 
 