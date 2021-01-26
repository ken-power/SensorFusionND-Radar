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
FMCW Waveform Design | Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp. | For given system requirements the calculated slope should be around 2e13 | PLANNED
Simulation Loop | Simulate Target movement and calculate the beat or mixed signal for every timestamp. |A beat signal should be generated such that once range FFT implemented, it gives the correct range i.e the initial position of target assigned with an error margin of +/- 10 meters. | PLANNED
Range FFT (1st FFT) | Implement the Range FFT on the Beat or Mixed Signal and plot the result. | A correct implementation should generate a peak at the correct range, i.e the initial position of target assigned with an error margin of +/- 10 meters. | PLANNED
2D CFAR | Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map. | The 2D CFAR processing should be able to suppress the noise and separate the target signal. The output should match the image shared in walkthrough. | PLANNED
2D CFAR | Create a CFAR README File | In a README file, write brief explanations for the following: <ul><li>Implementation steps for the 2D CFAR process.</li> <li>Selection of Training, Guard cells and offset. </li> <li>Steps taken to suppress the non-thresholded cells at the edges.</li> </ul> | [IN PROGRESS](#PRoject-Report)


# Project Report

## Implementation steps for the 2D CFAR process

## Selection of Training, Guard cells and offset

## Steps taken to suppress the non-thresholded cells at the edges


 