from math import *
import numpy as np
import matplotlib.pyplot as plt

from routines import *
from dust_opacity import *
from Primary import *
from Grid import *
from LambdaOp import *
from Detector import *
from EulerRotation import *

# 2

wave = 1.6
k=1          
r_grid = np.linspace(8.251*const.AU, 30*const.AU, 25*k)
phi_grid = np.linspace(-pi, pi, 36*k)
z_grid = np.linspace(-7*const.AU, 7*const.AU, 11*k)

grid = Grid(r_grid, phi_grid, z_grid)

for node in grid.nodelist:
   node.temperature = 1500.0

ori = EulerRotation(0,19,0)

det = Detector(200, 1.57, 1*50)
IMG = det.get_image(grid, 1.6, ori)

im = plt.imshow(IMG, interpolation='none', origin='lower')
bar = plt.colorbar(im, aspect=10, shrink=0.5)

plt.show()













