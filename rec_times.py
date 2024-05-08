"""
rec_times.py
============

Calculates the recurrence times for the chosen recurrence region.
Upper: I0 < 0
Lower: I0 > 0

Usage
-----
This script takes on two parameters
        python rec_times.py eps I0

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 08/05/2024
"""
import numpy as np
import sys
from functions import *
import os
# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)
# --- User input --- #
eps = float(sys.argv[1]) # The perturbation
I0 = float(sys.argv[2]) # The initial action.
# If I0 > 0, the script calculates the recurrence times for the LOWER region.
# If I0 < 0, the script calculates the recurrence times for the UPPER region.
# --- Parameters of the simulation --- #
N = int(1e9) # Maximum iteration time
# Variables to format the datafile
eN = int(np.log10(N))
bN = int(N/10**eN)
Irec = -np.sign(I0)
# Initial angle
theta0 = 0.5
# Calculates the recurrence times
rec_times = recurrence_times(theta0, I0, eps, Irec, N)
# Saves the data
df = "%s/recurrence_times_eps=%.5f_I0=%.5f_Irec=%i_N=%ie%i.dat" % (path, eps, I0, Irec,  bN, eN)
np.savetxt(df, rec_times)