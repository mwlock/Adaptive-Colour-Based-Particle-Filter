<h1 align="center">Object Tracking using a Color-Based Particle Filter with an Adaptive Target Color Distribution</h1>

|   |   |
|---|---|
| [<img src="/images/miguel_report.png" width="100%">](https://www.overleaf.com/read/qbrvkbphsbhq) <br> Miguel's Report| [<img src="/images/matthew_report.png" width="100%">](https://www.overleaf.com/read/ydqcxkjznwsd) <br> Matthew's Report |

**Authors** 
> [Miguel Garcia Naude](https://github.com/migsdigs), [Matthew Lock](https://github.com/matthew-william-lock)

## Abstract

Several methods exist to track objects. The particle filter has proven to be useful in non-linear and non-Gaussian tracking problems. The use of color features for object is advantageous as it can provide a framework that is rotation and scale invariant and robust with respect to partial-occlusion of an object and noise. However, color feature based algorithms are sensitive in particular to illumination and appearance changes of the tracked object.  A review and implementation of an existing object tracking algorithm using a color-based particle filter is conducted. This algorithm attempts to solve the problem of illumination and object appearance changes by adapting the color distribution of the target during gradual changes. Tracking is achieved by comparing the color distribution of a target to that of hypothetical states (particles). Color distributions are represented using histograms in the HSV color space, in an attempt to reduce the algorithms sensitivity to illumination changes. It is shown that color-based tracking is sufficient for use in ideal environments where the target has a distinct color distribution, and is robust to non-rigidity of objects, as well as partial occlusion. However, it appears as though methods to address slowly a adapting target are insufficient in cases where noise and occlusion is slowly introduced and is ineffective in cases of sudden illumination and appearance changes.

## Excerpt of Results

|   |   |
|---|---|
| ![car_Static](/simulation_results/gifs/car_static.gif) <br> Static target|  ![car_dynamic](/simulation_results/gifs/car_dynamic.gif) <br> Adaptive target |
|![freethrow](/simulation_results/gifs/freethrow.gif) <br> Simple object tracking|  ![surveillance](/simulation_results/gifs/surveillance.gif) <br> Adaptive lighting confitions|

An excerpt of the results shown above allude to the advantages of dynamic colour targets for objects with slowly changing color distributions. Other results highlight the filters general ability to track, as well as highlighting the shortfalls of tracking under extremely dynamic lighting conditions.

## Running the code

Further documentation regarding the execution of the simulations produced for this project can be found in the ```examples``` folder of this repository.
