"""
survival_probability.py
=======================

Calculates the survival probability given the limits of the survival region,
the number of initial conditions, and the maximum iteration time.

Usage
-----
        python survival_probability.py I_esc n_ic Nmax

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 2024-05-14
"""
import numpy as np
from functions import *
from datetime import datetime
import sys
import os

# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)

# --- User input --- #
I_esc = float(sys.argv[1]) # Limits of the survival region
n_ic = int(float(sys.argv[2])) # Number of initial conditions
Nmax = int(float(sys.argv[3])) # The maximum iteration time
# --- Parameters of the simulation --- #
eps = 1e-3 # The perturbation
# The lower and upper bound for theta
theta_ini = 0
theta_end = 1
# The lower and upper bound for I
dI = 1e-10
I_ini = -dI
I_end = dI
# Arrays with randomly chosen values for theta and I
theta = theta_ini + (theta_end - theta_ini)*np.random.rand(n_ic)
I = I_ini + (I_end - I_ini)*np.random.rand(n_ic)
# Variables to format the datafile
exponent = int(np.log10(n_ic))
base = int(n_ic/10**exponent)
# Calculates the escape times
esc_times = escape_time(theta, I, eps, I_esc, Nmax)
# Calculates the survival probability
sp = survival_prob(esc_times, Nmax)
# Saves the data
time = np.arange(1, Nmax + 1)
data = np.zeros((time.shape[0], 2))
data[:, 0] = time
data[:, 1] = sp
df = "%s/survival_probability_eps=%.5f_Iesc=%.10f_nic=%ie%i.dat" % (path, eps, I_esc, base, exponent)
np.savetxt(df, data)