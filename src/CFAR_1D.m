%% Implement 1D CFAR using lagging cells on the given noise and target scenario.

% Close and delete all currently open figures
close all;

%%

% Data_points - number of samples for which  we want to run the CFAR
Ns = 1000;

% Generate random noise for the same amount of samples
s=abs(randn(Ns,1));

%Targets location. 
% Assign random amplitudes [8, 9, 4, 11] to the noise samples at 
% bin 100, 200, 300, and 700 to make them appear as mock Targets
s([100 ,200, 300, 700])=[8 9 4 11];

%plot the signal to see the noise output
plot(s);

%% Apply CFAR to detect the targets by filtering the noise.

% 1. Define the necessary variables:

% 1a. Training Cells
T = 12;

% 1b. Guard Cells 
G = 4;

% Offset : Adding room above noise threshold for desired SNR 
% Here, we are working on linear values, so we multiply the offset
% by the threshold value.
offset=5;

% Vector to hold threshold values 
threshold_cfar = [];

%Vector to hold final signal after thresholding
signal_cfar = [];


%% Apply the 1D CIFAR

% 2. Slide window across the signal length
% step through all the cells from one end to another of the signal vector
% and make sure that we have the right spacing from the end

for i = 1:(Ns-(G+T+1))     

    % 2. - 5. Determine the noise threshold by measuring it within the training cells
    
    % for each step, add the noise within all the training cells
    noise_level = sum(s(i:i+T-1));
    
    % To determine the threshold take the average of the summed noise
    % and multiply it by the offset
    threshold = (noise_level/T) * offset;
    threshold_cfar = [threshold_cfar, {threshold}];

    % Pick the cell udner test, which is T+G cells away from the first
    % training cell and measure the signal level
    signal = s(i+T+G);
    
    
    % 6. Measuring the signal within the CUT

    % If the signal level at Cell Under Test is below the threshold
    % then assign it a 0 value
    if(signal < threshold)
        signal = 0;
    end

    % 8. Filter the signal above the threshold
    signal_cfar = [signal_cfar, {signal}];
end




% plot the filtered signal
plot (cell2mat(signal_cfar),'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')
