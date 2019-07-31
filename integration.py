from math import *
import numpy as np
import units_cgs as units
import phys_const_cgs as const
import matplotlib.pyplot as plt
from DiskStructure import *
from dust_opacity import *

def test():


   xx = np.linspace(8.3, 100, 300)

   dd = [struct.dust_density((x*const.AU, 0, 0)) for x in xx]
   TT = [struct.approx_temperature(x*const.AU) for x in xx]

   plt.plot(xx, dd, '.-')
   plt.show()

wavelength = 1.6
OP = DustOpacity()
kappa = OP.inter(wavelength)
print('kappa=', kappa, 'cm^2/g')

struct = DiskStructure()
dtau = 0.1
x1 = 8.3*const.AU
xx=[x1]
for i in range(50):
   print('x=', x1/const.AU, 'AU')
   rho1 = struct.dust_density((x1, 0, 0))
   dx = dtau/kappa/rho1
   print('dx=', dx/const.AU, 'AU')
   rho2 = struct.dust_density((x1+dx, 0, 0))
   dx2 = 2*dtau/kappa/(rho1+rho2)
   print('dx2=', dx2/const.AU, 'AU')
   print('-----')
   x1 += dx2
   xx.append(x1)


plt.plot([x/const.AU for x in xx])
plt.show()






