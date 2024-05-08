"""
escape_time_and_basin.py
========================

Calculates the escape time and escape basin for a grid of initial conditions.

Usage
-----
This script takes on two parameters.

        python escape_time_and_basin.py I_esc grid

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 08/05/2024
"""

import numpy as np
from joblib import Parallel, delayed
from functions import *
import sys
import os
# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)

# --- User input --- #
I_esc = float(sys.argv[1]) # The limits of the survival region
grid = int(sys.argv[2]) # The grid size
eps = 1e-3 # The perturbation
Nmax = int(1e6) # The maximum iteration time
# Lower and upper bound for theta
theta_ini = 0
theta_end = 1
# An array of uniformly distributed values
theta = np.linspace(theta_ini, theta_end, grid)
# Lower and upper bound for I
I_ini = -I_esc
I_end = I_esc
# An array of uniformly distributed values
I = np.linspace(I_ini, I_end, grid)
# Creates a grid of initial conditions
theta, I = np.meshgrid(theta, I)
# Calculates the escape time and escape basin
escape = np.array(Parallel(n_jobs=-1)(delayed(escape_time_and_basin)(theta[i, j], I[i, j], eps, I_esc, Nmax) for i in range(grid) for j in range(grid)))
esc_time = escape[:, 0].reshape((grid, grid))
esc_side = escape[:, 1].reshape((grid, grid))
# Writes the data into the file df
df = "%s/escape_times_eps=%.5f_Iesc=%.3f_grid=%i.dat" % (path, eps, I_esc, grid)
with open(df, "w") as df:
    for i in range(grid):
        for j in range(grid):
            df.write("%.16f %.16f %i %i\n" % (theta[i, j], I[i, j], esc_time[i, j], esc_side[i, j]))
        df.write("\n")