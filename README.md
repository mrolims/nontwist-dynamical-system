# nontwist-dynamical-system

Code repository accompanying the publication entitled "Ratchet current and scaling invariance in a nontwist mapping".

This project contains the code to generate and plot the data from all figures.

## Requirements

The required packages are listed in ``` requirements.txt ```. To install them please execute ``` pip install -r requirements.txt ```.

## Figure 1

To generate the fixed points, phase space, and plot the figure, run all cells within the heading named ``` Fig. 1 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 2

To generate the data of Figure 2, run ``` python lyapunov_vs_eps.py arg1 arg2 ```. Here, ``` arg1 ``` is the sample size of ``` eps ``` and ``` arg2 ``` is the total iteration time. Use ``` arg1 = 1000 ``` and ``` arg2 = 1e8 ```. To plot the figure, run all cells within the heading named ``` Fig. 2 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 3

To generate Figure 3, run all cells within the heading named ``` Fig. 3 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 4

To generate the data used in Figure 4, run ``` python exe_survival_probability.py arg ```. Here, ```arg ``` corresponds to each Figure, i.e., ``` arg = a ``` generates the data of Figure 4(a), ``` arg = b ``` generates the data of Figure 4(b), and ``` arg = c ``` generates the data of figure 4(c). To create the plot and perform the power law fitting, run all cells within the heading named ``` Fig. 4 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 5

To generate the data of Figure 5, run ``` python escape_time_and_basin.py arg ```. Here, ``` arg ``` corresponds to the limits of the survival region, ``` I_esc ```. You need to run it four times, for the values of ``` I_esc ``` mentioned in the figure caption. To create the plot and obtain the values in Table II, run all cells within the heading named ``` Fig. 5 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 6

To generate the data of Figure 6, run ``` python escape_side_ratio_vs_Iesc.py arg```. Here, ``` arg ``` corresponds to the perturbation ``` eps ```. You need to run it three times using the values of ``` eps ``` mentioned in the figure caption. To create the plot, run all cells within the heading named ``` Fig. 6 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 7

### Figure 7(a)

To generate the data of Figure 7(a), compile the Fortran program ``` I_averages.f90 ``` with either [``` gfortran ```](https://fortran-lang.org/learn/os_setup/install_gfortran/) or [``` ifx ```](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.8j44s3). For example, ``` ifx functions.f90 I_averages.f90 -o Iav.x ```, and run it as ``` ./ Iav.x arg1 arg2 arg3 ```. Here, ``` arg1 ```, ``` arg1 ```, and ``` arg1 ``` corresponds to the perturbation ``` eps ```, the number of initial conditions ``` n_ic ```, and the maximum iteration time ``` N ```. For Figure 7(a), you need to run it using the values of ``` eps ``` mentioned in the figure with ``` n_ic = 1e6 ``` and ``` N = 1e7 ```.

### Figures 7(b) and 7(c)

To generate the data of Figures 7(b) and 7(c), run ``` python rec_times arg1 arg2 ```. Here, ``` arg1 ``` and ``` arg2 ``` corresponds to the perturbation ``` eps ``` and the initial action ``` I0 ```, respectively. You need to run it with the values of ``` eps ``` mentioned in the figure with ``` I0 = 1e-10 ``` and ``` I0 = -1e-10 ``` to calculate the recurrence times of the lower and upper regions, respectively.


### Creating the plot

To create the plot and generate the data of Table III, run all cells within the heading named ``` Fig. 7 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 8

To generate the data of Figure 8, run the auxiliary script ``` exe_Irms.py 1e4 1e8 ```. It will compile the ``` I_averages.f90 ``` program with ``` ifx ``` and execute it for several values of ``` eps ``` (save them, you will need them for Figure 9). To create the plot, run all cells within the heading named ``` Fig. 8 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Figure 9

To generate the data of Figure 9, run the auxiliary script ``` exe_Irms.py 1e4 1e8 ```. It will compile the ``` I_averages.f90 ``` program with ``` ifx ``` and execute it for several values of ``` eps ```. To create the plot and calculate the critical exponents, run all cells within the heading named ``` Fig. 9 ``` in the ``` Plots.ipynb ``` Jupyter notebook.

## Contact

[matheusrolim95@gmail.com](mailto:matheusrolim95@gmail.com)
