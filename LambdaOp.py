import os
from math import *
import numpy as np
from routines import *
from DiskStructure import *
from dust_opacity import *
from bisect import bisect
import phys_const_cgs as const
import units_cgs as units
import scipy.sparse as sps
from scipy.integrate import simps
from radiation_transfer import DustRadTrans
from blackbody import B_lambda
from Primary import *
from Grid import *


class DirectionGrid:
   def __init__(self, n):
      self.grid = [np.array([0,0,1]), np.array([0,0,-1])]
      dphi = pi/n
      thetas = list(np.linspace(0, pi, n))
      del(thetas[-1])
      del(thetas[0])
      for theta in thetas:
         rc = sin(theta)
         L = 2*pi*rc
         m = floor(L/dphi)
         phis = list(np.linspace(-pi, pi, m))
         del(phis[-1])
         for phi in phis:
            x = sin(theta) * cos(phi)
            y = sin(theta) * sin(phi)
            z = cos(theta)
            self.grid.append(np.array([x,y,z]))


class LambdaOp:
   def __init__(self, wavelength, grid):
      self.struct = DiskStructureBench(True)
      self.opacity = DustOpacity()
      self.wave = wavelength
      self.kappa = self.opacity.inter(self.wave) + self.opacity.inter_sca(self.wave)
      print('kappa=', self.kappa)
      self.grid = grid
      self.dir_grid = DirectionGrid(7)
      self.dtau = 0.1
      self.tau_n = 20
      self.dx_max = 1*const.AU
      self.compute_density_integrals(100)

   def propagate_ray(self, x0, direct, dtau, n):
      tau_max = self.tau_n * self.dtau
      ray = Ray(x0, direct)
      dx_sum = 0
      tau = 0.0
      points = [x0]
      taus = [tau]
      #for i in range(n-1):
      while tau < tau_max:
         rho1 = self.struct.dust_density(ray.point_at(dx_sum))
         #if rho1==0: break
         dx1 = dtau/self.kappa/rho1
         #if isinf(dx1): print('inifnity')
         dx1 = min(dx1, self.dx_max)
         p_next = ray.point_at(dx_sum + dx1)
         p_cyl = CartCoord(p_next[0], p_next[1], p_next[2]).to_cyl()
         if not self.grid.is_point_inside(p_cyl): break
         rho2 = self.struct.dust_density(p_next)
         dx2 = 2*dtau/self.kappa/(rho1+rho2)
         dx2 = min(dx2, self.dx_max)
         dx_sum += dx2
         rho2 = self.struct.dust_density(ray.point_at(dx_sum))
         drho = dx2 * 0.5 * (rho1 + rho2)
         dtau = self.kappa * drho
         tau += dtau
         points.append( ray.point_at(dx_sum) )
         taus.append(tau)
      return points, taus

   def get_local_coords(self, point):
      p_cyl = CartCoord(point[0], point[1], point[2]).to_cyl()
      NN = self.grid.find_NN(p_cyl)
      neigh = self.grid.get_neighbors(NN)

      if neigh.r2!=None:
         hr = neigh.r2.cyl_coords.r - NN.cyl_coords.r
      else:
         hr = NN.cyl_coords.r - neigh.r1.cyl_coords.r

      x = (p_cyl.r - NN.cyl_coords.r)/hr

      h_phi_1 = neigh.phi2.cyl_coords.phi - NN.cyl_coords.phi
      h_phi_2 = NN.cyl_coords.phi - neigh.phi1.cyl_coords.phi
      h_phi = min(h_phi_1, h_phi_2)

      y = (p_cyl.phi - NN.cyl_coords.phi)/h_phi

      if neigh.z2!=None:
         hz = neigh.z2.cyl_coords.z - NN.cyl_coords.z
      else:
         hz = NN.cyl_coords.z - neigh.z1.cyl_coords.z

      z = (p_cyl.z - NN.cyl_coords.z)/hz
      if not -1.0<x<1.0: print('Warning: x=', x)
      if not -1.0<y<1.0: print('Warning: y=', y)
      if not -1.0<z<1.0: print('Warning: z=', z)
      return x,y,z, NN, neigh

   def get_source_at(self, point):
      x,y,z,NN,neigh = self.get_local_coords(point)
      S0 = 1 - x**2 - y**2 - z**2
      Sxp = 0.5*(x**2 + x)
      Sxm = 0.5*(x**2 - x)
      Syp = 0.5*(y**2 + y)
      Sym = 0.5*(y**2 - y)
      Szp = 0.5*(z**2 + z)
      Szm = 0.5*(z**2 - z)
      nodes = [NN, neigh.r2, neigh.r1, neigh.phi2, neigh.phi1, neigh.z2, neigh.z1]
      SF = []
      for node in nodes:
         if node==None:
            SF.append(0.0)
         else:
            SF.append( node.source )
      coeffs = [S0, Sxp, Sxm, Syp, Sym, Szp, Szm]
      S = sum([ z[0]*z[1] for z in zip(coeffs, SF) ])
      return S

   def transfer_row(self, node):
      #print('node:', node)
      #print('neighbors:', self.grid.get_neighbors(node))
      row = {}
      for direct in self.dir_grid.grid:
         x0 = node.cyl_coords.to_cart()
         ray_points, taus = self.propagate_ray(np.array([x0.x, x0.y, x0.z]), direct, self.dtau, self.tau_n)
         print('ray length:', len(ray_points))
         tau = 0.0
         i = 0
         for point in ray_points:
            x,y,z, NN, neigh = self.get_local_coords(point)

            S0 = 1 - x**2 - y**2 - z**2
            Sxp = 0.5*(x**2 + x)
            Sxm = 0.5*(x**2 - x)
            Syp = 0.5*(y**2 + y)
            Sym = 0.5*(y**2 - y)
            Szp = 0.5*(z**2 + z)
            Szm = 0.5*(z**2 - z)
            tau += taus[i]
            q = exp(-tau) * taus[i]
            i += 1

            set_or_add(row, NN.nid, q*S0)
            if neigh.r2 != None:
               set_or_add(row, neigh.r2.nid, q*Sxp)
            if neigh.r1 != None:
               set_or_add(row, neigh.r1.nid, q*Sxm)

            set_or_add(row, neigh.phi2.nid, q*Syp)
            set_or_add(row, neigh.phi1.nid, q*Sym)

            if neigh.z2 != None:
               set_or_add(row, neigh.z2.nid, q*Szp)
            if neigh.z1 != None:
               set_or_add(row, neigh.z1.nid, q*Szm)
      return row

   def compute_density_integrals(self, Npoints):
      if self.struct.symmetric:
         for i in range(self.grid.grid.shape[0]):
            for k in range(self.grid.grid.shape[2]):
               node = self.grid.grid[i,0,k]
               D = self.struct.dust_density_cyl(node.cyl_coords.r, node.cyl_coords.z, node.cyl_coords.phi)
               I = self.struct.density_integral(node.cyl_coords, Npoints)
               for j in range(self.grid.grid.shape[1]):
                  self.grid.grid[i,j,k].density_integral = I
                  self.grid.grid[i,j,k].tau = I*self.kappa
                  self.grid.grid[i,j,k].density = D
      else:
         for i in range(self.grid.grid.shape[0]):
            for j in range(self.grid.grid.shape[1]):
               for k in range(self.grid.grid.shape[2]):
                  node = self.grid.grid[i,j,k]
                  D = self.struct.dust_density_cyl(node.cyl_coords.r, node.cyl_coords.z, node.cyl_coords.phi)
                  I = self.struct.density_integral(node.cyl_coords, Npoints)
                  node.density_integral = I
                  node.tau = I*self.kappa
                  node.density = D


   def transfer_matrix(self):
      dg = len(self.dir_grid.grid)
      coo_rows = []
      coo_cols = []
      coo_data = []
      L = len(self.grid.nodelist)
      i=0
      for node in self.grid.nodelist:
         row = self.transfer_row(node)
         for el in row.items():
            coo_rows.append(node.nid)
            coo_cols.append(el[0])
            coo_data.append(el[1]/dg)
         nz = len(row)
         i += 1
         Tmax = max(row.values())
         Tsum = sum(row.values())
         node.Tmax = Tmax
         node.Tsum = Tsum
      COO = sps.coo_matrix( (coo_data, (coo_rows, coo_cols)), shape=(L,L) )
      K = COO.tocsr()   
      return K

   def precompute_matrices(grid, folder):
      ww = list(np.logspace(log10(0.5), log10(50), 15))
      ww.append(1.6)
      for w in ww:
         print('computing', w)
         op = LambdaOp(w, grid)
         T = op.transfer_matrix()
         sps.save_npz(folder+'\\T_'+'%.4e_um'%w+'.npz', T)
         print(w, 'computed')































