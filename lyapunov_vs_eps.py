"""
lyapunov_vs_eps.py
==================

Calculates the largest Lyapunov exponent as a function of eps.

Usage
-----
This script takes on two parameters

        python lyapunov_vs_eps.py n N

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 2024-05-14
"""
import numpy as np
from functions import *
from joblib import Parallel, delayed
import sys
import os
# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)

# --- User input --- #
n = int(sys.argv[1]) # Number of eps samples
N = int(float(sys.argv[2])) # Maximum iteration time
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