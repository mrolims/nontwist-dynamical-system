"""
lyapunov_vs_eps.py
==================

Calculates the largest Lyapunov exponent as a function of eps.

Usage
-----
        python lyapunov_vs_eps.py

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 08/05/2024
"""
import numpy as np
from functions import *
from joblib import Parallel, delayed
import os
# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)

# --- Parameters of the simulation --- #
n = 1000 # Number of eps samples
N = int(1e3) # Maximum iteration time
# Initial condition
psi = 0.5
I = 1e-10
# Variables to format the datafile
exponent = int(np.log10(N))
base = int(N/10**exponent)
# Array with the values of eps distributed in a log scale
eps = np.logspace(-6, -2, n)
# Calculates the largest Lyapunov exponent
lypnv = np.array(Parallel(n_jobs=-1)(delayed(lyapunov)(psi, I, eps[j], N) for j in range(n)))
lypnv = lypnv[:, 0]
# Saves the data
data = np.zeros((n, 2))
data[:, 0] = eps
data[:, 1] = lypnv
df = "%s/lyapunov_vs_eps_N=%ie%i_n=%i.dat" % (path, base, exponent, n)
np.savetxt(df, data)