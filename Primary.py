from math import *
import numpy as np
import units_cgs as units
import phys_const_cgs as const
import kepler as orb


class Primary:
   Rstar = 85*const.Rsol
   Lstar = 17500*const.Lsol
   Tstar = 7250

   ecc = 0.22
   P = 506
   sma = 1.2*const.AU
   i = np.radians(19)
   LAN = np.radians(0)

   def position(self, phase):
      "phase: 0 - periastron, 0.5 - apoastron"
      op = self
      E = orb.ecc_anomaly(op.ecc, phase * op.P, 0, op.P)
      theta = orb.true_anomaly(E, op.ecc)
      r_star = orb.radius(op.sma, op.ecc, theta)
      X_star,Y_star = r_star*cos(theta), r_star*sin(theta)
      return (X_star,Y_star, 0)

if __name__=='__main__':
   pri = Primary()
   print('%.2e %.2e %.2e'%pri.position(0))
   print('%.2e %.2e %.2e'%pri.position(0.5))



