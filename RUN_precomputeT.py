from Grid import *
from LambdaOp import *

grid = Grid.prep_grid_bench()

if bool(0):
   op = LambdaOp(0.55, grid)
   grid.write_VTK('out\\_density.vtk')
   print('density VTK saved')
   exit()

folder = r'D:\astro\Master\data\T_bench'
LambdaOp.precompute_matrices(grid, folder)



