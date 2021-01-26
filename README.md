# Radar Target Generation and Detection

The goal of this project is to use MATLAB to implement a Radar target generation and detection system.

![](images/project-layout.png)

This involves a number of steps, as shown in the diagram above:
* Configure the FMCW waveform based on the system requirements.
* Define the range and velocity of target and simulate its displacement.
* For the same simulation loop process the transmit and receive signal to determine the beat signal
* Perform Range FFT on the received signal to determine the Range
* Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.


## Project Specification

Function | Criteria | Specification | Status
:--- | :--- | :--- | :--- 
**FMCW Waveform Design** | Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp. | For given system requirements the calculated slope should be around 2e13 | PLANNED
**Simulation Loop** | Simulate Target movement and calculate the beat or mixed signal for every timestamp. |A beat signal should be generated such that once range FFT implemented, it gives the correct range i.e the initial position of target assigned with an error margin of +/- 10 meters. | PLANNED
**Range FFT (1st FFT)** | Implement the Range FFT on the Beat or Mixed Signal and plot the result. | A correct implementation should generate a peak at the correct range, i.e the initial position of target assigned with an error margin of +/- 10 meters. | PLANNED
**2D CFAR** | Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map. | The 2D CFAR processing should be able to suppress the noise and separate the target signal. The output should match the image shared in walkthrough. | PLANNED
**2D CFAR** | Create a CFAR README File | In a README file, write brief explanations for the following: <ul><li>Implementation steps for the 2D CFAR process.</li> <li>Selection of Training, Guard cells and offset. </li> <li>Steps taken to suppress the non-thresholded cells at the edges.</li> </ul> | [IN PROGRESS](#PRoject-Report)


# Project Report

## Source Code
The main source code is in the MATLAB file [radar_target_generation_and_detection.m](src/radar_target_generation_and_detection.m).


## System Requirements

### Radar System Requirements
The sensor fusion design for different driving scenarios requires different system configurations from a Radar. In this project, we design a Radar based on the following system requirements:

Requirement | Value
:--- | ---:
Frequency | 77 GHz
Range resolution | 1 m
Max range | 200 m
Max velocity | 70 m/s
Velocity resolution | 3 m/s

Max Range and Range Resolution will be considered here for waveform design.

* The sweep bandwidth can be determined according to the range resolution and the sweep slope is calculated using both sweep bandwidth and sweep time.


_Bandwidth(Bsweep)=speedoflight / (2∗rangeResolution)_

* The sweep time can be computed based on the time needed for the signal to travel the unambiguous maximum range. In general, for an FMCW radar system, the sweep time should be at least 5 to 6 times the round trip time. This example uses a factor of 5.5.

_Tchirp = 5.5 * 2 * Rmax / c_

Giving the slope of the chirp signal

_Slope=Bandwidth / Tchirp_


### Initial Range and velocity of the Target
We need to provide the initial range and velocity of the target. Range cannot exceed the max value of 200m and velocity can be any value in the range of -70 to + 70 m/s.


## Implementation steps for the 2D CFAR process

## Selection of Training, Guard cells and offset

## Steps taken to suppress the non-thresholded cells at the edges


 