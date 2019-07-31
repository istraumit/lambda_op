from math import *
import numpy as np

class EulerRotation:
   def __init__(self, alpha_deg, beta_deg, gamma_deg):
      "LAN, incl, arg. of  periastron"
      alpha = np.radians(alpha_deg)
      beta = np.radians(beta_deg)
      gamma = np.radians(gamma_deg)
      c1 = cos(alpha)
      c2 = cos(beta)
      c3 = cos(gamma)
      s1 = sin(alpha)
      s2 = sin(beta)
      s3 = sin(gamma)
      M = np.zeros((3,3))
      M[0, 0] = c1 * c3 - c2 * s1 * s3
      M[0, 1] = -c1 * s3 - c2 * c3 * s1
      M[0, 2] = s1 * s2
      M[1, 0] = c3 * s1 + c1 * c2 * s3
      M[1, 1] = c1 * c2 * c3 - s1 * s3
      M[1, 2] = -c1 * s2
      M[2, 0] = s2 * s3
      M[2, 1] = c3 * s2
      M[2, 2] = c2
      self.R = M

   def rotate(self, vector):
      assert(type(vector)==np.ndarray)
      assert(vector.shape==(3,))
      return self.R.dot( vector )

if __name__=='__main__':
   eu = EulerRotation(45, 20, 0)
   V = np.array([0,0,1])
   print(eu.R)
   print(eu.rotate(V))










