% Calculate the range of four targets with respective measured beat frequencies [0 MHz, 1.1 MHz, 13 MHz, 24 MHz].
% Given: The radar maximum range = 300m and the range resolution = 1m.

% Find the Bsweep of chirp for 1 m resolution
c = 3 * 10^8;                   % speed of light in meters/sec
delta_R = 1;                    % range resolution in meters
Bsweep = c/(2 * delta_R);     

% Calculate the chirp time based on the Radar's Max Range
range_max = 300;                % radar's max range
Ts = 5.5 * (range_max*2 / c);   % 5.5 times the trip time for maximum range

% define the frequency shifts 
beat_freq = [0 1.1e6 13e6 24e6];    % given beat frequencies for all four targets
calculated_range = c * Ts * beat_freq / (2 * Bsweep);

% Display the calculated range
disp(calculated_range);


% This outputs the following:
% 0   12.1000  143.0000  264.0000
% 
% indicating the first target is 0 meters away, the second target is 
% 12.1 meters away, the third is 143 meters away, and the fourth 
% is 264 meters away.