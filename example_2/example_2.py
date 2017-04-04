import numpy as np
import matplotlib.pyplot as plt
from beam import beam
from beam import npUtils
try:
    from beam import pyUtils
except:
    print("pyUtils not imported")
try:
    from beam import cyUtils
except:
    print("cyUtils not imported")
try:
    from beam import pcUtils
except:
    print("pcUtils not imported")

# Parameters
E  = 1e8
L  = 1.0
h0 = 2e-2
h1 = 1e-2
P0 = 1.0
P1 = 2.0
I0 = E*h0**3/12

# Creating the solver object
opts = {
    "boundary"  : ["cc", "ff"],
    "length"    : L,
    "nElem"     : 1000,
    "modulus"   : E
    }
sol = beam.FESolver(**opts)

# Constant loading, constant thickness
f1 = lambda x: P0/I0/24.0 * x**2 * (x**2 - 4.0*L*x + 6.0*L*L)
sol(utils = npUtils,
    load = lambda x: P0*np.ones_like(x),
    thickness = [h0])
g = sol.plotSolution(sty='b-')
plt.plot(sol.coor, f1(sol.coor), 'k:', label='Case 1')

plt.show()
