
from scipy.optimize import differential_evolution, minimize
import subprocess
import re
import numpy as np



def f_target(r1, r2, e, l, meshSizeFactor):
    exe = "./main"
    args = [str(r1) , str(r2), str(e), str(l), str(meshSizeFactor)]
    process = subprocess.Popen([exe] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    file = open("out.txt", "r")
    freq = file.readline()
    freq = float(freq)
    file.close()

    print("\n\n====Iteration Resume====")
    print(f"r1 = {r1}")
    print(f"r2 = {r2}")
    print(f"e = {e}")
    print(f"l = {l}")
    print(f"frequency = {freq}")
    return freq

r1  = 6e-3
r2  = 8e-3
e  = 20e-3
l   = 50e-3
meshSizeFactor = 0.6



bds = [(1e-3, 9e-3), (10e-3, 15e-3), (20e-3, 80e-3) ,(50e-3, 200e-3)]
target = 2420

f = lambda x: (f_target(x[0], x[1], x[2], x[3], meshSizeFactor) - target)**2 

#res = differential_evolution(f, bounds= bds, tol = 10-5)
#print(res)
res = minimize(f, x0 = [r1, r2, e, l], bounds = bds, tol = 1e-5, method = "Nelder-Mead")
print(res)