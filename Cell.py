from math import *
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg.dsolve as dsl
from routines import *
import phys_const_cgs as const

class Boundary:
   def set_flux(self, flux):
      self.incoming_flux = flux

   def __init__(self):
      self.incoming_flux = 0.0

   def __str__(self):
      return 'BND'

class Cell:
   def __init__(self, Nid, D, X, dX):
      assert(type(X)==CylCoord)
      assert(type(dX)==CylCoord)
      self.Nid = Nid
      self.D = D
      self.X = X
      self.dX = dX

   def min_dim(self):
      d2 = self.X.r * self.dX.phi
      return min(self.dX.r, d2, self.dX.z)

   def __str__(self):
      return str(self.Nid)

   def set_neighbors(self, neighbors):
      assert(type(neighbors)==Neighbors)
      self.neighbors = neighbors

   def get_row(self):
      c = const.c
      row = []
      Ah_sum = 0

      nei = self.neighbors.r1
      A = self.dX.phi*(self.X.r - 0.5*self.dX.r)*self.dX.z
      if type(nei)==Cell:
         h = cyl_dist(self.X, nei.X)
         Ah = A/h
         Ah_sum += Ah
         row.append( (nei.Nid, -c * self.D * Ah) )
      elif type(nei)==Boundary:
         h = self.dX.r
         Ah_sum += A/h

      nei = self.neighbors.r2
      A = self.dX.phi*(self.X.r + 0.5*self.dX.r)*self.dX.z
      if type(nei)==Cell:
         h = cyl_dist(self.X, nei.X)
         Ah = A/h
         Ah_sum += Ah
         row.append( (nei.Nid, -c * self.D * Ah) )
      elif type(nei)==Boundary:
         h = self.dX.r
         Ah_sum += A/h

      A = self.dX.r * self.dX.z
      for nei in [self.neighbors.phi1, self.neighbors.phi2]:
         if type(nei)==Cell:
            h = cyl_dist(self.X, nei.X)
            Ah = A/h
            Ah_sum += Ah
            row.append( (nei.Nid, -c * self.D * Ah) )
         elif type(nei)==Boundary:
            h = self.X.r * self.dX.phi
            Ah_sum += A/h

      A = self.X.r * self.dX.r * self.dX.phi
      for nei in [self.neighbors.z1, self.neighbors.z2]:
         if type(nei)==Cell:
            h = cyl_dist(self.X, nei.X)
            Ah = A/h
            Ah_sum += Ah
            row.append( (nei.Nid, -c * self.D * Ah) )
         elif type(nei)==Boundary:
            h = self.dX.z
            Ah_sum += A/h
      
      row.append( (self.Nid, c * self.D * Ah_sum) )
      return row

   def get_RHS(self):
      RHS = 0
      nei = self.neighbors.r1
      A = self.dX.phi*(self.X.r - 0.5*self.dX.r)*self.dX.z
      if type(nei)==Boundary:
         RHS += nei.incoming_flux * A

      nei = self.neighbors.r2
      A = self.dX.phi*(self.X.r + 0.5*self.dX.r)*self.dX.z
      if type(nei)==Boundary:
         RHS += nei.incoming_flux * A

      A = self.dX.r * self.dX.z
      for nei in [self.neighbors.phi1, self.neighbors.phi2]:
         if type(nei)==Boundary:
            RHS += nei.incoming_flux * A

      A = self.X.r * self.dX.r * self.dX.phi
      for nei in [self.neighbors.z1, self.neighbors.z2]:
         if type(nei)==Boundary:
            RHS += nei.incoming_flux * A

      return RHS



















   

