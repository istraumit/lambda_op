from math import *
from routines import *
import numpy as np

class Geometry:
   def intersect_plane(ray:Ray, plane_position:float, plane:int)->np.ndarray:
      "ray: object of Ray type"
      "plane_position: float"
      "plane: int (0-x, 1-y, 2-z)"
      if (ray.v[plane] == 0): return None
      t = (plane_position - ray.x0[plane]) / ray.v[plane]
      return ray.point_at(t)


   def intersect_cylinder(ray:Ray, cyl_radius:float, cyl_half_length:float)->tuple or None:

      def cyl_bounded(p, cyl_radius):
         return p[0]**2 + p[1]**2 < cyl_radius**2

      inter1 = Geometry.intersect_plane(ray, cyl_half_length, 2)
      inter2 = Geometry.intersect_plane(ray, -cyl_half_length, 2)
      if ray.v[0]**2 + ray.v[1]**2 == 0:
         if cyl_bounded(ray.x0, cyl_radius):
            return (inter1, inter2)
         else:
            return None

      A = (ray.v[0])**2 + (ray.v[1])**2
      B = 2 * (ray.x0[0] * ray.v[0] + ray.x0[1] * ray.v[1])
      C = (ray.x0[0])**2 + (ray.x0[1])**2 - (cyl_radius)**2

      D = QuadraticEquation.discriminant(A, B, C)
      if (D < 0): return None
      (t1,t2) = QuadraticEquation.solve(A, B, D)
      y1 = ray.x0[2] + t1 * ray.v[2]
      y2 = ray.x0[2] + t2 * ray.v[2]


      if (abs(y1) > cyl_half_length):
          if (abs(y2) > cyl_half_length):
              if (y1 * y2 < 0):
                  return (inter1, inter2)
              else: return None
          else:
              if (cyl_bounded(inter1, cyl_radius)):
                  return (ray.point_at(t2), inter1)
              else:
                  return (ray.point_at(t2), inter2)
      else:
          if (abs(y2) > cyl_half_length):
              if (cyl_bounded(inter2, cyl_radius)):
                  return (ray.point_at(t1), inter2)
              else:
                  return (ray.point_at(t1), inter1)
          else:
              return (ray.point_at(t1), ray.point_at(t2))

if __name__=='__main__':
   ray = Ray( np.array([-1.1,0,0]), np.array( [cos(pi/4), 0, sin(pi/4)] ) )
   x = Geometry.intersect_cylinder(ray, cyl_radius=1.0, cyl_half_length=1.0)
   print(x)





















