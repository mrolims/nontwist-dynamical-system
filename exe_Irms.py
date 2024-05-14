"""
exe_Irms.py
===========

Auxiliary script to run the I_averages.f90 program

Usage
-----
This scripts takes on two parameters

        python exe_Irms n_ic N

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 2024-05-14
"""
import os
import sys
import numpy as np

n_ic = int(float(sys.argv[1]))
N = int(float(sys.argv[2]))


eps = np.logspace(-4, -2, 20)
I_ini = 1e-10
I0 = np.ones_like(eps, dtype=np.float64)*I_ini

f90_file = "I_averages.f90"
exe_file = "Irms.x"

print("Removing the executable file to avoid old file execution...")
comm = "rm -rf %s" % exe_file
print("$", comm)
os.system(comm)
print("Succesfully removed %s.\nCompiling files..." % exe_file)
comm = "ifx functions.f90 %s -o %s" % (f90_file, exe_file)
print("$", comm)
os.system(comm)

if os.path.isfile(exe_file):
    print("Compilation succeeded. Executing the program...")
else:
    print("Compilation failed. Aborting execution.")
    import sys
    sys.exit()

for i in range(eps.shape[0]):
    comm = "./%s %.16f %i %i" % (exe_file, eps[i], n_ic, N)
    print(comm)
    os.system(comm)
