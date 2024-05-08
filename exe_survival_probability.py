"""
exe_survival_probability.py
===========================

Auxiliary script to run the survival_probability.py script.

Usage
-----
This scripts takes on one parameter

        python survival_probability.py fig

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 08/05/2024
"""
import os
import numpy as np
import sys

# To obtain the data needed for Fig. 4(a), run with fig = a
# To obtain the data needed for Fig. 4(b), run with fig = b
# To obtain the data needed for Fig. 4(c), run with fig = c
fig = sys.argv[1]
# Number of initial conditions
n_ic = int(1e6)

if fig == "a":
    I_escs = np.arange(0.01, 0.045, 0.005)
elif fig == "b":
    I_escs = np.logspace(np.log10(1e-3), np.log10(4e-2), 100)
elif fig == "c":
    I_escs = np.arange(0.002, 0.011, 0.001)

for Iesc in I_escs:
    comm = "python survival_probability.py %.10f %i" % (Iesc, n_ic)
    print("$", comm)
    os.system(comm)
