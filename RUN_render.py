import phys_const_cgs as const
from Grid import *
from EulerRotation import *
from Detector import *
import matplotlib.pyplot as plt
import numpy as np

name = 'symm'

SOURCE = []
with open('out\\'+name+'_solution.dat') as f:
   for line in f:
      row = line.split()
      SOURCE.append(float(row[0]))

print('solution loaded:', len(SOURCE), 'nodes')

k=1          
r_grid = np.linspace(8.251*const.AU, 30*const.AU, 25*k)
phi_grid = np.linspace(-pi, pi, 36*k)
z_grid = np.linspace(-7*const.AU, 7*const.AU, 11*k)

grid = Grid(r_grid, phi_grid, z_grid)
for i in range(len(SOURCE)):
   grid.nodelist[i].source = SOURCE[i]

print('grid ready')

inclin = 19
LAN = -90
img_num = 0
#while LAN<=90:
for LAN in [-90, 90]:
   print('rendering LAN='+str(LAN)+'...')
   ori = EulerRotation(LAN, inclin, 0)

   det = Detector(200, 1.57, 0.5*50)
   IMG = det.get_image(grid, 1.6, ori)
   np.save('out\\image_'+str(LAN), IMG)

   im = plt.imshow(IMG, interpolation='none', origin='lower')
   bar = plt.colorbar(im, aspect=20, shrink=1.0)
   bar.set_label('Intensity [erg/cm^2/s/cm/sr]')
   plt.tight_layout()
   num = str(img_num).zfill(2)
   plt.savefig('out\\'+num+'_figure_'+str(LAN)+'.pdf')
   plt.savefig('out\\'+num+'_figure_'+str(LAN)+'.png')
   plt.clf()
   LAN += 10
   img_num += 1


























