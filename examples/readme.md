# Object tracking with colour histogram based Particle filter

This file is used to run the implemented particle filter for tracking. Inluded in this file are 4 example video sequences as described in the project paper. The examples can be run by opening the ``` particle_filter_tracking.m ``` file and following the instructions provided below.

## Running the code

1. Add data to the Matlab environment by right clicking on the ```data``` folder in the main repository directory and selecting ```Add to path > Selected Folder and Subfolders```
2. Add source code to the Matlab environment by right clicking on the ```src``` folder in the main repository directory and selecting ```Add to path > Selected Folder and Subfolders```
3. Switch between video sequences by changing the ```simulation``` parameters as described in 'particle_filter_tracking.m'
4. Adapt other parameters as you please (including turning on/off localisation, chaning process noise, etc). Instructions are provided in the 'particle_filter_tracking.m' file.

## Chaning the number of bins for colour histogram

The number of bins is set to {16x16x1} for the HSV colour space by default. This can be changed by navigating to the ```get_weighted_histogram.m``` file and changing the number of bins for each channel. The file is location in ```\src\matlab\colour_functions\get_weighted_histogram.m```