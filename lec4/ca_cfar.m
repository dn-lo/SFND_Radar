% Implement 1D CFAR using lagging cells on the given noise and target scenario.

%% Generation
% Close and delete all currently open figures
close all;

% Data_points
Ns = 1000;

% Generate random noise
s=abs(randn(Ns,1));

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s([100 ,200, 300, 700])=[8 9 4 11];

%plot the output

%% Apply CFAR to detect the targets by filtering the noise.

% 1. Define the following:
% 1a. Training Cells (half-amplitude t)
% 1b. Guard Cells (half-amplitude g)
T = 12; 
G = 2;

% Offset : Adding room above noise threshold for desired SNR 
offset = 3;

% Vector to hold threshold values 
threshold_cfar = zeros(Ns,1);

%Vector to hold final signal after thresholding
signal_cfar = zeros(Ns,1);

% 2. Slide window across the signal length
for i = 1:Ns   

    % 2. - 5. Determine the noise threshold by measuring it within the training cells
    ist = max(1, i - (G + T));                  % window starting point: must be greater than initial position
    iend = min(Ns, i + (G + T));                % window ending point: must be smaller than final position
    ids = [ist:(i-G-1) (i+G+1):iend];           % indices to use for noise calculation
    threshold_cfar(i) =  offset * mean(s(ids)); % compute mean then use offset for SNR
    
    % 6. Measuring the signal within the CUT
    % 8. Filter the signal above the threshold
    if (s(i) > threshold_cfar(i))
        signal_cfar(i) = s(i);
    else
        signal_cfar(i) = 0;
    end
end


% plot the filtered signal
% plot (signal_cfar,'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(threshold_cfar,'r--','LineWidth',2)
hold on, plot(signal_cfar,'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')