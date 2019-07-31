from math import *
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg.dsolve as dsl
import phys_const_cgs as const
from routines import *
from Cell import *
from DiskStructure import *
from Rosseland import *

class Grid:
   beta = [] # ratios of the Rosseland mean free path to the min dimension of the cell

   def __init__(self, nphi, r1, r2, nr, z1, z2, nz):
      rr = np.logspace(log10(r1), log10(r2), nr+1)
      #rr = np.linspace(r1, r2, nr+1)
      zz = np.linspace(z1, z2, nz+1)
      phis = np.linspace(-pi, pi, nphi+1)

      disk = DiskStructure()
      ross = Rosseland(np.linspace(10, 2000, 300))

      D = 1.0
      grid = np.empty( (nphi, nr, nz), dtype=Cell )
      nid=0
      for i in range(nphi):
         for j in range(nr):
            for k in range(nz):
               X = CylCoord(0.5*(rr[j] + rr[j+1]), 0.5*(phis[i]+phis[i+1]), 0.5*(zz[k]+zz[k+1]))
               dX = CylCoord(rr[j+1]-rr[j], phis[i+1]-phis[i], zz[k+1]-zz[k])
               T = disk.approx_temperature(X.r)
               rho = disk.dust_density(cyl2cart(X))
               kappa_R = ross.inter(T)
               assert(rho>0)
               assert(kappa_R>0)
               D = 1./(3 * rho * kappa_R)
               cell = Cell(nid, D, X, dX)
               self.beta.append(3*D/cell.min_dim())
               nid += 1
               grid[i,j,k] = cell

      for i in range(nphi):
         for j in range(nr):
            for k in range(nz):
               phi1 = grid[i-1, j, k]
               if i+1<nphi: phi2 = grid[i+1, j, k]
               else:        phi2 = grid[0, j, k]
               
               if j-1<0: r1 = Boundary()
               else: r1 = grid[i, j-1, k]

               if j+1>=nr: r2 = Boundary()
               else: r2 = grid[i, j+1, k]

               if k-1<0: z1 = Boundary()
               else: z1 = grid[i, j, k-1]

               if k+1>=nz: z2 = Boundary()
               else: z2 = grid[i, j, k+1]

               N = Neighbors(r1, r2, phi1, phi2, z1, z2)
               grid[i,j,k].set_neighbors(N)

      self.cells = grid

   def solve(self):
      print('assembling')
      L = len(self.cells.flat)

      coo_rows = []
      coo_cols = []
      coo_data = []
     
      RHS_ = []
      
      for cell in self.cells.flat:
         row = cell.get_row()
         for el in row:
            coo_rows.append(cell.Nid)
            coo_cols.append(el[0])
            coo_data.append(el[1])
         
         RHS_.append(cell.get_RHS())

      COO = sps.coo_matrix( (coo_data, (coo_rows, coo_cols)), shape=(L,L) )

      K = COO.tocsr()   
      RHS = np.asarray(RHS_)

      print('solving')
      U = dsl.spsolve(K, RHS, use_umfpack=False)
      print('solved')

      for cell in g.cells.flat:
         cell.U = U[cell.Nid]


import time

n = 30
nphi = 25
nr = 30
nz = 2

start = time.time()
g = Grid(nphi, 10*const.AU, 50*const.AU, nr, -1*const.AU, 1*const.AU, nz)

print('Ross_MFP / R_cell_min:', min(g.beta), max(g.beta))

exit()

for i in range(0, nphi//2):
   for j in range(nz):
      cell = g.cells[i,0,j]
      cell.neighbors.r1.set_flux(10.0)

end = time.time()
print('Grid prep time:', end-start, 's')


start = time.time()
g.solve()
end = time.time()
print('Solver time:', end-start, 's')


exit()

import matplotlib.pyplot as plt

IMG = np.zeros((nphi, nr))
for i in range(nphi):
   for j in range(nr):
      IMG[i,j] = g.cells[i,j,0].U


im = plt.imshow(IMG, interpolation='none', origin='lower')
bar = plt.colorbar(im, aspect=10, shrink=0.5)

plt.show()



IMG = np.zeros((nr, nz))
for i in range(nr):
   for j in range(nz):
      IMG[i,j] = g.cells[5,i,j].U


im = plt.imshow(IMG, interpolation='none', origin='lower')
bar = plt.colorbar(im, aspect=10, shrink=0.5)

plt.show()


















   

