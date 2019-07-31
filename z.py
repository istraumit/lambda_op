from math import *
import inspect
import numpy as np
import matplotlib.pyplot as plt

R_peak = 25
p1 = 2.0
p2 = 1.0

def sigma(R):
   if R<R_peak:
      return (R/R_peak)**p1
   else:
      return (R_peak/R)**p2


L = 2
H = 200
IMG = np.zeros((H,H))

for i in range(H):
   for j in range(H):
      x = i - H//2
      y = j - H//2
      rxy = sqrt(x**2 + y**2)
      phi = atan2(y,x)
      if rxy<1: rho = 0.0
      else: rho = 0.5*(sin(L*phi + 20*log10(rxy)) + 1) * sigma(rxy)
      if rho<0: rho = 0.0
      IMG[i,j] = rho


im = plt.imshow(IMG, interpolation='none', origin='lower', extent=[-H//2, H//2, -H//2, H//2])

bar = plt.colorbar(im, aspect=20, shrink=1)

plt.xlabel('X [AU]')
plt.ylabel('Y [AU]')
plt.show()

















      