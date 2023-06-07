# Guide through files
written by Arnaud Hallemans 2022-2023 at IMEC/NERF

## Assembly detector
Contains an automatic assembly detector based on literature mentioned in file. 
Input is the binned neuronal activity throughout learning. Outputs the number of estimated assemblies and which neurons form those assemblies and at which times those assemblies are active.

FastICA_25 is downloaded from https://github.com/aludnam/MATLAB/tree/master/FastICA_25

## Detection algorithm
Identifying neurons in 4 types.

First run extract_after_spikes.m with corresponding spiking files in the right directories, then let classigying_after_stim.m run for each patient (saving and loading files from extract_after_spike.m).

## Dorsal vs Ventral
up_down_regulation computes the difference in activity between late and early phase of learning and performs several statistical tests on them.

First run extract_after_spikes.m, then up_down_regulation.m.


## Correlation between height and spiking rate
Visualisation of the correlations between foot height and spiking rates.

First run extract_height_fr.m, then visualise_height_fr.m. tfix.m can be used to align the stimulations to the height data (have to include those files yourself)


## Sleep classifier
Classifies sleep based on EEG and EMG data into REM, NREM and wake.
