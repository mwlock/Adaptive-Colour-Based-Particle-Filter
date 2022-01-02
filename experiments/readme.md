# Results

This section documents the experiments that were performed and the best parameters found for each experiment.

## HSV Colour Histogram

| Parameter 	| Description                  	| Value 	|
|-----------	|------------------------------	|-------	|
| h_bins    	| Number of bins for H channel 	| 16    	|
| s_bins    	| Number of bins for H channeS 	| 16    	|
| v_bins    	| Number of bins for V channel 	| 1     	|

## Dynamically tracking a basketball with an ellipsoidal region

| Parameter                          	| Description                              	| Value                   	|
|------------------------------------	|------------------------------------------	|-------------------------	|
| Scale                              	| Image scale                              	| 0.3                     	|
| M                                  	| Number of particles                      	| 400                     	|
| R                                  	| Process covariance matrix                	| diag([20 20 5 5])*scale 	|
| mean_state_observation_prob_thresh 	| Observation probability threshold        	| 0.45                    	|
| alpha                              	| Contribution of the mean state histogram 	| 0.05                    	|

- We were able to successfully localize with 400 particles. Less than 300 results in a failed localisation.

## Dynamically tracking bee in a video game with an ellipsoidal region

| Parameter                          	| Description                              	| Value                   	|
|------------------------------------	|------------------------------------------	|-------------------------	|
| Scale                              	| Image scale                              	| 0.3                     	|
| M                                  	| Number of particles                      	| 200                     	|
| R                                  	| Process covariance matrix                	| diag([50 50 20 20])*scale |
| mean_state_observation_prob_thresh 	| Observation probability threshold        	| -                     	|
| alpha                              	| Contribution of the mean state histogram 	| -                     	|

- It was found that applying a dynamic target model to this experiment would result in the target distribution slowly ressembling the distribution of the background and ultimately failing to track the bee.
- Since the bee is quite small and hard to find, the pf was initialised to track the bee. This avoids needing an exorbitant amount of particles for localisation.

## Dynamically tracking car game with an ellipsoidal region

| Parameter                          	| Description                              	| Value                   	|
|------------------------------------	|------------------------------------------	|-------------------------	|
| Scale                              	| Image scale                              	| 0.3                     	|
| M                                  	| Number of particles                      	| 200                     	|
| R                                  	| Process covariance matrix                	| diag([50 50 20 20])*scale |
| mean_state_observation_prob_thresh 	| Observation probability threshold        	| -                     	|
| alpha                              	| Contribution of the mean state histogram 	| -                     	|

- It was found that applying a dynamic target model to this experiment would result in the target distribution slowly ressembling the distribution of the background. While tracking was still possible, the region size grew to an extent that dramatically increased the computational load.
- Since the car is quite small and hard to find, the pf was initialised to track the car. This avoids needing an exorbitant amount of particles for localisation.

## Surveillance tracking car game with an ellipsoidal region

| Parameter                          	| Description                              	| Value                   	|
|------------------------------------	|------------------------------------------	|-------------------------	|
| Scale                              	| Image scale                              	| 0.6                     	|
| M                                  	| Number of particles                      	| 200                     	|
| R                                  	| Process covariance matrix                	| diag([50 50 20 20])*scale |
| mean_state_observation_prob_thresh 	| Observation probability threshold        	| -                     	|
| alpha                              	| Contribution of the mean state histogram 	| -                     	|
