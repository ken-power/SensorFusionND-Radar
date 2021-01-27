clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 3e8; %speed of light
maxRange = 200;
rangeResolution = 1;
maxVelocity = 100;

%% User Defined Range and Velocity of target
% Define the target's initial position and velocity. 
% Note : Velocity remains contant
range = 110;
velocity = -20;


%% FMCW Waveform Generation

% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B_sweep = c/(2*rangeResolution);   % Bandwidth (B)
T_chirp = 5.5*2*maxRange/c;        % Chirp time
slope = B_sweep/T_chirp;           % Slope of the FMCW

disp("slope = ");
disp(slope);

% Operating carrier frequency of Radar 
fc = 77e9;             %carrier freq
                                                          
% The number of chirps in one sequence. Its ideal to have 2^ value for 
% the ease of running the FFT for Doppler Estimation. 
Nd = 128;      % # of doppler cells OR # of sent periods % number of chirps

% The number of samples on each chirp. 
Nr = 1024;     % for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t = linspace(0,Nd*T_chirp,Nr*Nd); % total time for samples


% Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t));    % transmitted signal
Rx = zeros(1,length(t));    % received signal
Mix = zeros(1,length(t));   % beat signal

% Similar vectors for range_covered and time delay.
r_t = zeros(1,length(t));
td = zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         

    % For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = range + velocity * t(i);
    td(i) = 2*r_t(i)/c;

    % For each time sample we need update the transmitted and
    % received signal. 
    Tx(i) = cos(2 * pi * (fc * t(i) + slope * (t(i)^2)/2));
    Rx(i)  = cos(2 * pi * (fc * (t(i) - td(i)) + slope * ((t(i)-td(i))^2)/2));
    
    % Now by mixing the Transmit and Receive generate the beat signal
    % This is done by element wise matrix multiplication of Transmit and
    % Receiver Signal
    Mix(i) = Tx(i) .* Rx(i);
    
end

%% RANGE MEASUREMENT - FFT Operation

% Implement the 1D FFT on the Mixed Signal

% Reshape the vector into Nr*Nd array. Nr and Nd here would also define 
% the size of Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr,Nd]);

% Run the FFT on the beat signal along the range bins dimension (Nr) and
% normalize.
sig_fft1 = fft(Mix,Nr);
sig_fft1 = sig_fft1./Nr;

% Take the absolute value of FFT output
sig_fft1 = abs(sig_fft1);

% Output of FFT is double sided signal, but we are interested in only one 
% side of the spectrum. Hence we throw out half of the samples, keeping
% one half of the signal.
single_side_sig_fft1 = sig_fft1(1:Nr/2);

% Plot the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% plot FFT output; 
% There should be a peak at the initial position of the target 
plot(single_side_sig_fft1); 
 
axis ([0 200 0 1]);


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map. Implement CFAR on the generated Range Doppler Map (RDM).


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix = reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

% use the surf function to plot the output of 2DFFT and to show axis 
% in both dimensions
doppler_axis = linspace(-100, 100, Nd);
range_axis = linspace(-200, 200, Nr/2)*((Nr/2)/400);
figure('Name', '2DFFT output - Range Doppler Map')
surf(doppler_axis, range_axis, RDM);


%% CFAR implementation

% Slide Window through the complete Range Doppler Map

% Determine the number of Training cells for each dimension. 
% Similarly, pick the number of guard cells.

% Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 8;

% Select the number of Guard Cells in both dimensions around the Cell under 
% test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
offset = 1.4;

% Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


% Design a loop such that it slides the CUT across range doppler map by
% giving margins at the edges for Training and Guard Cells.
% For every iteration sum the signal level within all the training
% cells. To sum convert the value from logarithmic to linear using db2pow
% function. Average the summed values for all of the training
% cells used. After averaging convert it back to logarithimic using pow2db.
% Further add the offset to it to determine the threshold. Next, compare 
% the signal under CUT with this threshold. If the CUT level > threshold
% assign it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR
RDM = RDM/max(max(RDM));

% The process above will generate a thresholded block, which is smaller 
% than the Range Doppler Map as the CUT cannot be located at the edges of
% matrix. Hence, few cells will not be thresholded. To keep the map size 
% same set those values to 0. 

   
% Slide the cell under test across the complete matrix. Make sure the 
% CUT has margin for Training and Guard cells from the edges.

for i = Tr+Gr+1:(Nr/2)-(Tr+Gr)
 
    for j = Td+Gd+1:(Nd)-(Td+Gd)
        
        % Create a vector to store noise_level for each iteration on 
        % training cells
        noise_level = zeros(1,1);
        
        % For every iteration sum the signal level within all the 
        % training cells. To sum convert the value from logarithmic to 
        % linear using db2pow function.
        
        % Step through each of bins and the surroundings of the CUT
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                
                % Exclude the Guard cells and CUT cells
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    % Convert db to power
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        
        % Calculate threshould from noise average then add the offset:

        % Average the summed values for all of the training cells used. 
        trainingCellsSummedValues = (2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1);
        
        % After averaging convert it back to logarithmic using pow2db        
        threshold = pow2db(noise_level / trainingCellsSummedValues);
        
        % Add the offset to determine the threshold
        threshold = threshold + offset;
        
        % Measure the signal in Cell Under Test(CUT) and compare against
        CUT = RDM(i,j);
        
        % Next, compare the signal under CUT against this threshold.
        % If the CUT level > threshold assign it a value of 1
        if (CUT < threshold)
            RDM(i,j) = 0;
        
        % else equate it to 0
        else 
            RDM(i,j) = 1;
        end
        
    end
end


% The process above will generate a thresholded block, which is smaller 
% than the Range Doppler Map as the CUTs cannot be located at the edges
% of the matrix due to the presence of Target and Guard cells. Hence, 
% those cells will not be thresholded.
%
% To keep the map size same as it was before CFAR, equate all the 
% non-thresholded cells to 0.
RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows
RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns 


% display the CFAR output using the Surf function like we did for Range
% Doppler Response output.
figure('Name','The output of the 2D CFAR process')
surf(doppler_axis,range_axis,RDM);
colorbar;



 
 