from math import *
import numpy as np
from routines import *
from bisect import bisect
from blackbody import B_lambda
import phys_const_cgs as const


class Node:
   def __init__(self, nid, grid_coords, cyl_coords):
      assert(type(cyl_coords)==CylCoord)
      self.nid = nid
      self.grid_coords = grid_coords
      self.cyl_coords = cyl_coords
   def __str__(self):
      return 'Node #'+str(self.nid)+':'+str(self.grid_coords)+str(self.cyl_coords)   
   def source_func(self, wave):
      return B_lambda(wave, self.temperature)


class Grid:
   def __init__(self, r_grid, phi_grid, z_grid):
      self.r_grid = r_grid
      self.phi_grid = phi_grid
      self.z_grid = z_grid
      grid = np.empty( (len(r_grid), len(phi_grid)-1, len(z_grid)), dtype=type(Node) )
      nodelist = []
      nid=0
      for i in range(len(r_grid)):
         for j in range(len(phi_grid)-1):
            for k in range(len(z_grid)):
               coords = CylCoord(r_grid[i], phi_grid[j], z_grid[k])
               node = Node(nid, (i,j,k), coords)
               grid[i,j,k] = node
               node.value = j
               nodelist.append(node)
               nid += 1

      self.grid = grid
      self.nodelist = nodelist

   def find_NN(self, p_cyl):
      assert(type(p_cyl) == CylCoord)
      r_grid = self.r_grid
      phi_grid = self.phi_grid
      z_grid = self.z_grid

      i = bisect(r_grid, p_cyl.r)
      if i==len(r_grid) or (i>0 and (abs(r_grid[i-1] - p_cyl.r)<abs(r_grid[i] - p_cyl.r))):
         i -= 1

      L_phi = self.grid.shape[1]
      j = bisect(phi_grid[:-1], p_cyl.phi)
      if j==0: pass
      elif j==L_phi: j -= 1
      elif (abs(phi_grid[j-1] - p_cyl.phi)<abs(phi_grid[j] - p_cyl.phi)): j -= 1 

      k = bisect(z_grid, p_cyl.z)
      if k==len(z_grid) or (k>0 and (abs(z_grid[k-1] - p_cyl.z)<abs(z_grid[k] - p_cyl.z))):
         k -= 1

      return self.grid[i,j,k]

   def get_neighbors(self, node):
      assert(type(node) == Node)
      r1, r2, z1, z2 = None, None, None, None
      (i,j,k) = node.grid_coords
      if i>0: r1 = self.grid[i-1, j, k]
      if i<len(self.r_grid)-1: r2 = self.grid[i+1, j, k]
      if k>0: z1 = self.grid[i, j, k-1]
      if k<len(self.z_grid)-1: z2 = self.grid[i, j, k+1]
      J = len(self.phi_grid)-1
      phi1 = self.grid[i, (j-1)%J, k]
      phi2 = self.grid[i, (j+1)%J, k]
      return Neighbors(r1, r2, phi1, phi2, z1, z2)   

   def is_point_inside(self, p_cyl):
      assert(type(p_cyl)==CylCoord)
      return self.r_grid[0]<p_cyl.r<self.r_grid[-1] and self.z_grid[0]<p_cyl.z<self.z_grid[-1]
      
   def write_VTK(self, fn):
      N = len(self.nodelist)
      with open(fn, 'w') as f:
         f.write('# vtk DataFile Version 3.0\n')
         f.write('title\n')
         f.write('ASCII\n')
         f.write('DATASET UNSTRUCTURED_GRID\n')
         f.write('POINTS '+str(N)+' float\n')
         for node in self.nodelist:
            cart_coords = cyl2cart(node.cyl_coords)
            f.write('%.4e %.4e %.4e\n'%cart_coords)
         
         g = self.grid
         cells = []
         for i in range(g.shape[0]-1):
            for j in range(g.shape[1]-1):
               for k in range(g.shape[2]-1):
                  cells.append( (g[i,j,k].nid, g[i+1, j, k].nid, g[i+1, j+1, k].nid, g[i, j+1, k].nid, \
                     g[i,j,k+1].nid, g[i+1, j, k+1].nid, g[i+1, j+1, k+1].nid, g[i, j+1, k+1].nid) )
                  
         for i in range(g.shape[0]-1):
            for k in range(g.shape[2]-1):
                  j1, j2 = 0, g.shape[1]-1
                  cells.append( (g[i,j1,k].nid, g[i+1, j1, k].nid, g[i+1, j2, k].nid, g[i, j2, k].nid, \
                     g[i,j1,k+1].nid, g[i+1, j1, k+1].nid, g[i+1, j2, k+1].nid, g[i, j2, k+1].nid) )
               
                  
                  
         f.write('CELLS '+str(len(cells))+' '+str(9*len(cells))+'\n')
         for cell in cells:
            s = str(cell).replace(',',' ')
            f.write('8 '+s[1:-1]+'\n')
         
         f.write('CELL_TYPES '+str(len(cells))+'\n')
         for i in range(len(cells)):
            f.write('12\n')
         
         f.write('POINT_DATA ' + str(N)+'\n')
         #fields = ['density', 'temperature', 'density_integral','tau','Tmax','Tsum']
         #fields.extend(['J_int', 'J_star_int', 'J_star', 'temperature_eq'])
         fields = []
         for dk in self.nodelist[0].__dict__:
            if type(self.nodelist[0].__dict__[dk])==np.float64:
               fields.append(dk)

         for field in fields:
            if hasattr(self.nodelist[0], field):
               f.write('SCALARS '+field+' float 1\n')
               f.write('LOOKUP_TABLE default\n')
               for i in range(N):
                  v = getattr(self.nodelist[i], field)
                  if v<1.e-35: v=0.0
                  f.write('%.8e'%( v )+'\n')

   def prep_grid_V390():
      k=1          
      r_grid = np.linspace(8.251*const.AU, 30*const.AU, 25*k)
      phi_grid = np.linspace(-pi, pi, 36*k)
      z_grid = np.linspace(-7*const.AU, 7*const.AU, 11*k)

      grid = Grid(r_grid, phi_grid, z_grid)
      return grid

   def prep_grid_bench():
      r_grid = np.logspace( log10(1*const.AU), log10(1000*const.AU), 50)
      #r_grid = np.linspace( 1*const.AU, 1000*const.AU, 50)
      phi_grid = np.linspace(-pi, pi, 18)
      z_grid = np.linspace(-15*const.AU, 15*const.AU, 11)

      grid = Grid(r_grid, phi_grid, z_grid)
      return grid











