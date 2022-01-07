# Object Tracking using a Color-Based Particle Filter with an Adaptive Target Color Distribution

## Abstract

Several methods exist to track objects. The particle filter has proven to be useful in non-linear and non-Gaussian tracking problems. The use of color features for object is advantageous as it can provide a framework that is rotation and scale invariant and robust with respect to partial-occlusion of an object and noise. However, color feature based algorithms are sensitive in particular to illumination and appearance changes of the tracked object.  A review and implementation of an existing object tracking algorithm using a color-based particle filter is conducted. This algorithm attempts to solve the problem of illumination and object appearance changes by adapting the color distribution of the target during gradual changes. Tracking is achieved by comparing the color distribution of a target to that of hypothetical states (particles). Color distributions are represented using histograms in the HSV color space, in an attempt to reduce the algorithms sensitivity to illumination changes. It is shown that color-based tracking is sufficient for use in ideal environments where the target has a distinct color distribution, and is robust to non-rigidity of objects, as well as partial occlusion. However, it appears as though methods to address slowly a adapting target are insufficient in cases where noise and occlusion is slowly introduced and is ineffective in cases of sudden illumination and appearance changes.

## Results

#### Simple target tracking

![freethrow](/simulation_results/gifs/freethrow.gif)
<!-- <img src="/simulation_results/gifs/freethrow.gif" width="40" height="40" /> -->

#### Adaptive tracking for motor vehicle

Static target            |  Adaptive target
:-------------------------:|:-------------------------:
![freethrow](/simulation_results/gifs/car_static.gif) |  ![freethrow](/simulation_results/gifs/simulation_2_adaptive.gif)