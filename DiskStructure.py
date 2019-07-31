from math import *
import numpy as np
import units_cgs as units
import phys_const_cgs as const
import matplotlib.pyplot as plt
from scipy.integrate import simps
from routines import *

class DiskParams:
   p1 = 2.0
   p2 = 1.0
   pT = 0.75
   sigma_max = 10.
   R_peak = 25*const.AU
   R_min = 2*const.AU
   R_max = 200*const.AU
   M_star = (0.79+1.86)*const.Msol
   R_sub = 1*const.AU # sublimation rim radius
   rim_shift = -0.4*const.AU*0 # sublimation rim shift towards apoastron due to long-term irradiation asymmetry
   T_sub = 1100.
   dust_to_gas = 0.01

class DiskStructure:
   params = DiskParams()
   RIM_CENTRE = (params.rim_shift, 0, 0)

   def __init__(self, symmetric):
      assert(type(symmetric)==bool)
      self.symmetric = symmetric

   def omega_K(self, R):
      return sqrt(const.G * self.params.M_star/R**3)

   def approx_temperature(self, R):
      return self.params.T_sub*(self.params.R_sub/R)**self.params.pT

   def sound_speed(self, R):
      T = self.approx_temperature(R)
      return sqrt( const.k*T/(2.3*const.mp) )

   def scale_height(self, R):
      c = self.sound_speed(R)
      om = self.omega_K(R)
      return c/om

   def sigma(self, R):
      params = self.params
      if R<params.R_peak:
         return params.sigma_max*(R/params.R_peak)**params.p1
      else:
         return params.sigma_max*(params.R_peak/R)**params.p2

   def azi_profile(self, phi):
      return 1.0 + sin(0.5*(phi + 0.5*pi))**4

   def spiral(self, phi, R):
      return 0.5*(sin(2*phi + 20*log10(R/const.AU)) + 1)

   def sigma_asymm(self, R, phi):
      return self.sigma(R) * self.spiral(phi, R) #self.azi_profile(phi)

   def gas_density(self, R, z, phi):
      H = self.scale_height(R)
      if self.symmetric:
         S = self.sigma(R)
      else:
         S = self.sigma_asymm(R, phi)
      rho_c = S/(sqrt(2*pi)*H)
      rho = rho_c * exp(-0.5*(z/H)**2)
      return rho

   def dust_density(self, p):
      (x,y,z) = p
      r, z, theta = cart_2_cyl(p)
      return self.dust_density_cyl(r, z, theta)

   def dust_density_cyl(self, r, z, phi):
      gas = self.gas_density(r, z, phi)
      d = sqrt(r**2 + z**2)
      if d>self.params.R_sub: return gas * self.params.dust_to_gas
      else: return 0.0

   def density_integral(self, p_cyl:CylCoord, Npoints):
      r, z, phi = p_cyl.r, p_cyl.z, p_cyl.phi
      r_sub = self.params.R_sub
      z0 = r_sub * z/r
      rr = np.linspace(r_sub, r, Npoints)
      zz = np.linspace(z0, z, Npoints)
      S = dist3d( (r_sub, z0, 0), (r, z, 0) )
      ss = np.linspace(0, S, Npoints)
      rho = [self.dust_density_cyl(z[0], z[1], phi) for z in zip(rr, zz)]
      I = simps(rho, ss)
      return I

class DiskStructureBench(DiskStructure):
   def dust_density_cyl(self, r, z, phi):
      kappa = 620 # cm^2/g
      tau = 1
      Rout = 1000*const.AU
      rd = 0.5*Rout
      zd = 0.25*rd
      f1 = rd/r
      h  = zd * (r/rd)**1.125
      f2 = exp( -(pi/4) * (z/h)**2 )
      rho0 = tau/(0.5 * kappa * Rout * log(Rout))
      rho = 4 * rho0 * f1 * f2
      return rho



if __name__=='__main__':

   struct = DiskStructureBench(True)
   H = 501
   IMG = np.zeros((H,H))
   extent = 100*const.AU
   pix_scale = extent/H

   for i in range(H):
      for j in range(H):
         x = pix_scale*(i - 0.5*H)
         y = pix_scale*(j - 0.5*H)

         IMG[i,j] = log10( 1 + struct.dust_density((y, 0, x)) )

   im = plt.imshow(IMG, interpolation='none', origin='lower')
   bar = plt.colorbar(im, aspect=20, shrink=1.0)

   plt.show()










