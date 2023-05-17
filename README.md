# Barrier Hamiltonian Monte Carlo (BHMC)

This is a Matlab implementation of Barrier Hamiltonian Monte Carlo (BHMC) based on the implementation of CRHMC [1] available at 

https://github.com/ConstrainedSampler/PolytopeSamplerMatlab

that comes with the submission "Unbiased constrained sampling with Self-Concordant Barrier Hamiltonian Monte Carlo".

[1] "Sampling with riemannian hamiltonian monte carlo in a constrained space", Kook et al., NEURIPS 2022

## Folder/File Structure

1. `Instances`: this folder includes constrained-based models from systems biology and LP test sets in `0raw` folder. Running `preparePolytope.m` prepares instances for each algorithm -- it preprocesses constrained-based models by a presolver in the CRHMC package, transforms them into full-dimensional instances, and then round the full-dimensional ones by the maximum volume ellipsoid (MVE) algorithm.
2. `BHMC-test`: `testpaper_bhmc.m` and `testpaper_rhmc.m` run BHMC and CRHMC respectively on the instances and stores result mat files in `bhmc_test` and `rhmc_test` subfolders.
3. `CHRR-test`: this folder contains CHRR files downloaded from Bounciness/Volume-and-Sampling repository. It is needed for some auxiliary functions.
5. `PolytopeSamplerMatlab-develop`: this folder contains the implementation of BHMC and CRHMC.
There are two strategies to choose the step size. 
The one that target an acceptance ratio of 0.5. The functions in the subfolder `code` implementing this strategy are `sample_modif.m` and `sample_GL.m`.
The second one is based on the module in the orignal code. The functions implementing this strategy are `sample.m` and `sample_GL_v2.m`.
6. `Benchmark`: `main_run_benchmark.m` runs BHMC and CRHMC on the examples of target distributions described in the paper. The results are stored mat files in subfolders. To plot the boxplots use the script `plot_bias_boxplot.m`.
7. `testpaper_plot.m`: based on sampling results from the bio and LP dataset, it draws plots for sampling time and mixing rate versus dimension and the number of nonzeros.
It also displays the sampling time per effective sample.


## Illustration of the bias
Run `main_run_benchmark.m` in the folder `Benchmark` The results are stored mat files in the subfolder `bias_analysis_local` that you should create if it does not exist. To plot the boxplots use the script `plot_bias_boxplot.m` in the same folder.

## Illustration on real datasets
### Data Preparation
You should first create two empty subfolders `1chrr` and `2cdhr` under `Instances` folder. Then simply run `preparePolytope.m` to prepare rounded instances for CHRR and CDHR.

### Algorithms on Bio and LP Instances
All algorithms are limited in their implementation to one core for fair comparison.

Run `testpaper_rhmc.m`and `testpaper_bhmc.m` in the folder `BHMC-test`
Plot the result using the function `testpaper_plot.m`


