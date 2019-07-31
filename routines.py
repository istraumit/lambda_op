from math import *
import numpy as np

class CylCoord:
   def __init__(self, r, phi, z):
      self.r = r
      self.phi = phi
      self.z = z
   def __str__(self):
      return '(%.2e, %.2e, %.2e)'%(self.r, self.phi, self.z)
   def to_cart(self):
      x = self.r * cos(self.phi)
      y = self.r * sin(self.phi)
      z = self.z
      return CartCoord(x,y,z)


class CartCoord:
   def __init__(self, x, y, z):
      self.x = x
      self.y = y
      self.z = z
   def __str__(self):
      return '(%.2e, %.2e, %.2e)'%(self.x, self.y, self.z)
   def to_cyl(self):
      r = sqrt(self.x**2 + self.y**2)
      phi = atan2(self.y, self.x)
      return CylCoord(r, phi, self.z)
   def asarray(self):
      return np.array([self.x, self.y, self.z])

class Ray:
   def __init__(self, x0, v):
      assert(type(x0)==np.ndarray)
      assert(type(v)==np.ndarray)
      assert( abs(np.linalg.norm(v)-1)<1.e-10 )
      self.x0 = x0
      self.v = v
   def point_at(self, dx):
      return self.x0 + dx*self.v

def cart_2_cyl(p):
   (x,y,z) = p
   r = sqrt(x**2 + y**2)
   theta = atan2(y,x)
   return r, z, theta

def dist3d(p1, p2):
   return sqrt( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 )

def cyl2cart(p):
   assert(type(p)==CylCoord)
   x = p.r * cos(p.phi)
   y = p.r * sin(p.phi)
   z = p.z
   return (x,y,z)

def cyl_dist(p1, p2):
   assert(type(p1)==CylCoord)
   assert(type(p1)==CylCoord)
   c1 = cyl2cart(p1)
   c2 = cyl2cart(p2)
   return dist3d(c1, c2)

# set or add element to dictionary
def set_or_add(D, i, v):
   if i in D: D[i] += v
   else: D[i] = v

class Neighbors:
   def __init__(self, r1, r2, phi1, phi2, z1, z2):
      self.r1 = r1
      self.r2 = r2
      self.phi1 = phi1
      self.phi2 = phi2
      self.z1 = z1
      self.z2 = z2

   def __str__(self):
      s = []
      s.append(str(self.r1))
      s.append(str(self.r2))
      s.append(str(self.phi1))
      s.append(str(self.phi2))
      s.append(str(self.z1))
      s.append(str(self.z2))
      return ','.join(s)

class QuadraticEquation:
   def discriminant(a,b,c):
      return b * b - 4 * a * c

   def solve(a,b,D):
      sq = sqrt(D)
      root1 = (-b + sq) / (2 * a)
      root2 = (-b - sq) / (2 * a)
      return (root1, root2)














