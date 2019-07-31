from math import *
import numpy as np

class CylCoord:
   def __init__(self, r, phi, z):
      self.r = r
      self.phi = phi
      self.z = z
 
   def __str__(self):
      return '(%.2e, %.2e, %.2e)'%(self.r, self.phi, self.z)

def cyl2cart(p):
   x = p.r * cos(p.phi)
   y = p.r * sin(p.phi)
   z = p.z
   return (x,y,z)

class Node:
   def __init__(self, nid, grid_coords, cyl_coords):
      self.nid = nid
      self.grid_coords = grid_coords
      self.cyl_coords = cyl_coords

 
class Grid:
   def __init__(self, r_grid, phi_grid, z_grid):
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
         f.write('SCALARS density float 1\n')
         f.write('LOOKUP_TABLE default\n')
         for i in range(N):
            f.write(str(np.random.randn())+'\n')
            #f.write(str(self.nodelist[i].value)+'\n')
         
k=1          
r_grid = np.logspace(1, 2, 25*k)
phi_grid = np.linspace(-pi, pi, 36*k)
z_grid = np.linspace(-20, 20, 11*k)

g = Grid(r_grid, phi_grid, z_grid)
g.write_VTK('vtk_test.vtu')




