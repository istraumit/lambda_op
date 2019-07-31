from math import *
import numpy as np
import phys_const_cgs as const
import units_cgs as units
from routines import *
from LambdaOp import *
from Geometry import *

class Detector:
   def __init__(self, size_pix, distance_kpc, FOV_mas):
      self.size_pix = size_pix
      self.distance_kpc = distance_kpc
      self.FOV_mas = FOV_mas

   def get_image(self, grid, wavelength_um, orient):

      def coincide(p1, p2):
         return dist3d(p1, p2)<1.0e3

      def order_points(inter_outer_):
         if inter_outer_[0][2] < inter_outer_[1][2]:
            inter_outer = inter_outer_
         else:
            inter_outer = (inter_outer_[1], inter_outer_[0])
         return inter_outer


      R_inner = grid.r_grid[0]
      R_outer = grid.r_grid[-1]
      halfLen = 0.5*(grid.z_grid[-1]-grid.z_grid[0])

      H = self.size_pix
      center = H/2.0
      pixel_scale = self.FOV_mas / H
      print('Pixel_scale:', pixel_scale)
      pixel_AU = pixel_scale * self.distance_kpc * const.AU

      op = LambdaOp(wavelength_um, grid)
      direct = orient.rotate( np.array([0,0,-1]) )
      dtau = 0.1
      ntau = 20
      tau_grid = []
      tau = 0
      for i in range(ntau):
         tau_grid.append(tau)
         tau += dtau

      IMG = np.zeros((H,H))
      for i in range(H):
         for j in range(H):
            IMG[i,j] = 0.0
            sky_X = (i - center)*pixel_AU
            sky_Y = (j - center)*pixel_AU
            X_local = orient.rotate(np.array([sky_X, sky_Y, 0]) )
            ray = Ray(X_local, direct)
            inter_outer_ = Geometry.intersect_cylinder(ray, R_outer, halfLen)
            inter_inner_ = Geometry.intersect_cylinder(ray, R_inner, halfLen)
            if inter_outer_==None: continue

            inter_outer = order_points(inter_outer_)


            if inter_inner_ == None:
               entry = inter_outer[1]
            else:
               inter_inner = order_points(inter_inner_)
               if coincide(inter_outer[0], inter_inner[0]) and coincide(inter_outer[1], inter_inner[1]):
                  continue
               else:
                  if inter_inner[0][2]>inter_outer[0][2]:
                     entry = inter_inner[0]
                     if inter_inner[1][2]<inter_outer[1][2]:
                        entry = inter_outer[1]
                  else:
                     entry = inter_outer[1]
            
            points = op.propagate_ray(entry, direct, dtau, ntau)
            SF = [op.get_source_at(point) for point in points]
            tau_grid2 = tau_grid[:len(SF)]
            integrand = [z[0]*exp(-z[1])*dtau for z in zip(SF, tau_grid2)]
            I = simps(integrand, tau_grid2)
            IMG[i,j] = I
      return IMG








      




