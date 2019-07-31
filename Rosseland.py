from math import *
import numpy as np
import matplotlib.pyplot as plt
import phys_const_cgs as co
import units_cgs as units
from dust_opacity import DustOpacity
from scipy.integrate import simps

class Rosseland:

   OP = DustOpacity()
   wave_grid = np.logspace( log10(0.1), log10(1.e4), 1000 )

   def dB_lambda_dT(lam, T):
      Q = co.h * co.c / lam / co.k / T
      if Q>100:
         Z = exp(-Q)
      else:
         Z = exp(Q)/(exp(Q) - 1)**2
      return (2 * co.h**2 * co.c**3 / lam**6 / co.k) * Z / T**2 

   def get_dB_dT_norm(lams, T):
      yy = [Rosseland.dB_lambda_dT(lam*units.um, T) for lam in lams]
      yy_max = max(yy)
      yy = [y/yy_max for y in yy]
      return yy

   def get_k_R(self, T):
      lams = self.wave_grid
      dB_dT = Rosseland.get_dB_dT_norm(lams, T)
      k_tot = self.OP.inter(lams)

      I_1 = simps([z[0]/z[1] for z in zip(dB_dT, k_tot)], lams)
      I_2 = simps(dB_dT, lams)
      k_R = I_2/I_1
      return k_R

   def inter(self, T):
      return np.interp(T, self.T_grid, self.k_R)

   def __init__(self, T_grid):
      self.T_grid = T_grid
      self.k_R = [self.get_k_R(T) for T in T_grid]

if __name__=='__main__':
   T_grid = np.linspace(10, 1500, 100)
   ROSS = Rosseland(T_grid)
   TT = np.linspace(0, 2000, 150)
   RR = ROSS.inter(TT)

   plt.plot(T_grid, ROSS.k_R, linewidth=2)
   plt.plot(TT, RR)
   plt.xlabel('Temperature [K]')
   plt.ylabel('Rosseland mean opacity [cm^2/g]')
   plt.show()










