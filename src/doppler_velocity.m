% Calculate the velocity in m/s of four targets with the following doppler
% frequency shifts: 
% 3KZ, -4.5KHz, 11KHz, -3KHz].
%
% Given: Radar's operating frequency is 77 GHz.

c = 3 * 10^8;       % speed of ligth
frequency = 77e9;   % radar operating frequency in Hz

% Calculate the wavelength
wavelength = c/frequency;

% Define the Doppler shifts in Hz
doppler_shift = [3e3 -4.5e3 11e3 -3e3];

% Calculate the velocity of the targets
% fd = 2 * Vr / lambda

Vr = doppler_shift * wavelength / 2;

disp(Vr);

% Output =  5.8442   -8.7662   21.4286   -5.8442
