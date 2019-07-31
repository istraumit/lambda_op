import os
from math import *
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg.dsolve as dsl

from routines import *
from dust_opacity import *
from Primary import *
from Grid import *
from LambdaOp import *
from Detector import *
from EulerRotation import *

wave = 1.6
opacity = DustOpacity()
epsilon = opacity.epsilon(wave)
print('Photon destruction prob.:', epsilon)

name = 'symm'

folder = r'D:\astro\Master\data\T_'+name
npzs = [f for f in os.listdir(folder) if f.endswith('.npz')]
npzs_sorted = sorted(npzs, key=lambda x:float(x.split('_')[1]) )
wave_grid_um = [float(x.split('_')[1]) for x in npzs_sorted]
wave_grid = [x*units.um for x in wave_grid_um]
#for i,v in enumerate(wave_grid_um): print(i,v)

print('Lambda-matrices')
for i,v in enumerate(npzs_sorted):
   print(i, v)

TT = []
for npz in npzs_sorted:
   TT.append(sps.load_npz(folder+'\\'+npz))


k=1          
r_grid = np.linspace(8.251*const.AU, 30*const.AU, 25*k)
phi_grid = np.linspace(-pi, pi, 36*k)
z_grid = np.linspace(-7*const.AU, 7*const.AU, 11*k)

grid = Grid(r_grid, phi_grid, z_grid)
op = LambdaOp(wave, grid)
sample_node = grid.grid[5, 18, 5]
#print(sample_node.cyl_coords)

T = TT[4]
U = sps.csr_matrix( np.diag([1 for x in range(len(grid.nodelist))]) )
Teps = U - (1-epsilon)*T

primary = Primary()
pri_coord = np.asarray( primary.position(0.5) )
RT = DustRadTrans(85*const.Rsol, 17000*const.Lsol)
for node in grid.nodelist:
   d = dist3d(cyl2cart(node.cyl_coords), pri_coord)
   node.J_star_int = RT.get_J_star_integral(d, node.density_integral)
   node.J_star = RT.get_J_star(d, node.density_integral, wave*units.um)
   node.temperature = RT.get_T_from_J(node.J_star_int)
   node.temperature_eq = node.temperature

   node_coord = node.cyl_coords.to_cart().asarray()
   pri_dir = node_coord - pri_coord
   norm = np.linalg.norm(pri_dir)
   pri_dir = (1/norm)*pri_dir
   node.primary_direction = pri_dir

print('zero iteration complete')

kappas = op.opacity.inter(wave_grid_um)

it = 0
print(it, ':', '%.1f'%sample_node.temperature)
while it<=10:
   JJ = []
   for i in range(len(TT)):
      wave = wave_grid[i]
      T = TT[i]
      S = [node.source_func(wave) for node in grid.nodelist]
      J = T.dot(S)
      JJ.append(J)

   for node in grid.nodelist:
      i = node.nid
      Js = [z[0][i]*z[1] for z in zip(JJ, kappas)]
      J_int = simps(Js, wave_grid)
      node.J_int = J_int + node.J_star_int
      node.temperature = RT.get_T_from_J(node.J_int)

   it += 1
   print(it, ':', '%.1f'%sample_node.temperature)

for node in grid.nodelist:
   node.RHS = (1-epsilon)*node.J_star + epsilon*node.source_func(wave*units.um)

RHS = np.asarray([node.RHS for node in grid.nodelist])

print('solving')
SOLU = dsl.spsolve(Teps, RHS, use_umfpack=False)
print('solved')

for i in range(len(SOLU)):
   grid.nodelist[i].source = SOLU[i]

grid.write_VTK('out\\'+name+'_model.vtk')
print('model saved')

with open('out\\'+name+'_solution.dat', 'w') as f:
   for node in grid.nodelist:
      f.write('%.6e %.6e %.6e %.6e %.6e'%(node.source, node.J_star,\
         node.primary_direction[0], node.primary_direction[1], node.primary_direction[2]))
      f.write('\n')
print('solution saved')










