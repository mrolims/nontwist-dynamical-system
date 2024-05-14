"""
escape_side_ratio_vs_Iesc.py
============================

Calculates the fraction of particles that escapes from either exit.
Also calculates some statistical measures regarding the escape times.

Usage
-----
This script takes on one parameters.

        python escape_side_ratio_vs_Iesc.py eps

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 2024-05-14
"""

import numpy as np
from functions import *
from joblib import Parallel, delayed
from collections import Counter
import sys
import os
# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)
# --- User input --- #
eps = float(sys.argv[1]) # The perturbation
#
if eps == 1e-3:
    I_esc_end = 0.045
elif eps == 2e-3:
    I_esc_end = 0.063
elif eps == 3e-3:
    I_esc_end = 0.074
# --- Parameters of the simulation --- #
Nmax = int(1e6) # Maximum iteration time
n_ic = int(1e6) # The number of initial conditions
I_esc_ini = 0.0001 # Lower bound for the I_esc
n_I_esc = 100 # Sample size of I_esc
I_esc = np.linspace(I_esc_ini, I_esc_end, n_I_esc) # Array with the values of I_esc uniformly distributed
theta_ini = 0 # Lower bound for theta
theta_end = 1 # Upper bound for theta
dI = 1e-10 # Lower and upper bound for the action
I_ini = -dI
I_end = dI
# Creates two array for the angle and the action with randomly chosen values
theta = theta_ini + (theta_end - theta_ini)*np.random.rand(n_ic)
I = I_ini + (I_end - I_ini)*np.random.rand(n_ic)
# Variables to format the datafile
exponent = int(np.log10(n_ic))
base = int(n_ic/10**exponent)
edI = int(np.log10(dI))
bdI = int(dI/10**edI)
# Files where the data will be stored
df = "%s/escape_side_ratio_eps=%.5f_dI=%ie%i_nic=%ie%i_nIesc=%i.dat" % (path, eps, bdI, edI, base, exponent, n_I_esc)
df_sr = open(df, "w")
df = "%s/escape_times_statistics_eps=%.5f_dI=%ie%i_nic=%ie%i_nIesc=%i.dat" % (path, eps, bdI, edI, base, exponent, n_I_esc)
df_ets = open(df, "w")
# Changes I_esc
for j in range(n_I_esc):
    escape = np.array(Parallel(n_jobs=-1)(delayed(escape_time_and_basin)(theta[i], I[i], eps, I_esc[j], Nmax) for i in range(n_ic)))
    esc_time = escape[:, 0]
    esc_side = escape[:, 1]
    #
    count = Counter(esc_side)
    esc_top = count[1]/n_ic
    island = count[0]/n_ic
    esc_bottom = count[-1]/n_ic
    df_sr.write("%.16f %.16f %.16f %.16f %.16f\n" % (I_esc[j], esc_top, island, esc_bottom, esc_top + island + esc_bottom))
    #
    bottom_index = np.where(esc_side == -1)
    top_index = np.where(esc_side == 1)
    esc_time_top = esc_time[top_index]
    esc_time_bottom = esc_time[bottom_index]
    measures = np.zeros(4)
    top_measures = np.zeros(4)
    bottom_measures = np.zeros(4)
    #
    measures[0] = np.mean(esc_time)
    measures[1] = np.var(esc_time)
    measures[2] = np.std(esc_time)
    measures[3] = np.max(esc_time) - np.min(esc_time)
    top_measures[0] = np.mean(esc_time_top)
    top_measures[1] = np.var(esc_time_top)
    top_measures[2] = np.std(esc_time_top)
    top_measures[3] = np.max(esc_time_top) - np.min(esc_time_top)
    bottom_measures[0] = np.mean(esc_time_bottom)
    bottom_measures[1] = np.var(esc_time_bottom)
    bottom_measures[2] = np.std(esc_time_bottom)
    bottom_measures[3] = np.max(esc_time_bottom) - np.min(esc_time_bottom)
    #
    df_ets.write("%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n" % (I_esc[j], measures[0], measures[1], measures[2], measures[3], top_measures[0], top_measures[1], top_measures[2], top_measures[3], bottom_measures[0], bottom_measures[1], bottom_measures[2], bottom_measures[3]))