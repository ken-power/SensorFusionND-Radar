%Operating frequency (Hz)
fc = 77.0e9;

%Transmitted power (W)
Ps = 3e-3;

%Antenna Gain (linear)
G =  10000;

%Minimum Detectable Power
Pe = 1e-10;

%RCS of a car
RCS = 100;

%Speed of light
c = 3*10^8;

%Calculate the wavelength
wavelength = c / fc;
disp("wavelength = " + wavelength);

%Measure the Maximum Range a Radar can see. 
R = ((Ps * G^2 * wavelength^2 * RCS) / (Pe * (4 * 3.14)^3))^(-4)

disp("range = " + R);
