% 2-D Transform
% The 2-D Fourier transform is useful for processing 2-D signals and other 2-D data such as images.
% Create and plot 2-D data with repeated blocks.

P = peaks(20);
X = repmat(P,[5 10]);
imagesc(X)

% Compute the 2-D Fourier transform of the data.  

signal_fft = fft2(X);

% Shift the zero-frequency component to the center of the output
signal_fft = fftshift(signal_fft);

% take the absolute value
signal_fft = abs(signal_fft);

% plot the resulting 100-by-200 matrix, which is the same size as X.
imagesc(signal_fft);
