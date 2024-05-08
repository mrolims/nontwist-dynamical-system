"""
Irms_sat_vs_eps.py
============================

Extracts the saturation value for the I_rms given the datafile containing I_rms x n.

Usage
-----
        python Irms_sat_vs_eps.py

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 08/05/2024
"""
import numpy as np
from functions import *
import pandas as pd
import os
# The path variable defines the path to where the data will be stored.
path = "Data"
if not os.path.isdir(path):
    os.system("mkdir %s" % path)
# --- Parameters of the simulation --- #
N = int(1e8) # The maximum iteration time
n_ic = int(1e4) # The number of initial conditions
eps = np.logspace(-4, -2, 20) # The values of eps
I0 = 1e-10 # The initial acition
I_sat = np.zeros((20, 2)) # An array to store the eps and I_sat
# Variables to format the datafile
eN = int(np.log10(N))
bN = int(N/10**eN)
en_ic = int(np.log10(n_ic))
bn_ic = int(n_ic/10**en_ic)
eI0 = int(np.floor(np.log10(I0)))
bI0 = int(I0/10**eI0)

for i in range(len(eps)):

    df = "%s/Irms_eps=%.5f_I0=%ie%i_nic=%ie%i_N=%ie%i.dat" % (path, eps[i], bI0, eI0, bn_ic, en_ic, bN, eN)
    print("Extraindo dados de %s..." % df)
    df = pd.read_csv(df, header=None, delim_whitespace=True)
    time = np.array(df[0])
    Irms = np.array(df[1])
    I_sat[i, 1] = Irms[-1]

df = "%s/Isat_vs_eps_I0=%ie%i_nic=%ie%i_N=%ie%i.dat" % (path,bI0, eI0, bn_ic, en_ic, bN, eN)
I_sat[:, 0] = eps
np.savetxt(df, I_sat)