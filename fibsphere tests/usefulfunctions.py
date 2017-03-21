from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass
try: range = xrange
except NameError: pass

import numpy as np
import datetime
import time

# if a function will be used in multiple parts of the project, define it here.

# == GLOBAL CONSTANTS ==
c = 9.715611713408621e-12 # light speed in kpc/s
microG_TO_T = 1e-10
GeV_TO_kg = 1.78266184e-27
m_TO_kpc = 3.24077929e-20
kgkpc2_PER_s2_TO_GeV = 5.9427943e48

def dot(u, v):
    #u and v must be in cartesian coordinates
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

def mag(V):
    return np.sqrt(V[0]**2 + V[1]**2 + V[2]**2)

def theta_cylindrical(X, Y): # this theta goes from 0 to 2pi
    if X == 0.0:
        if Y > 0.0:
            return np.pi/2
        if Y < 0.0:
            return 3.0*np.pi/2
        else:
            return 0
    else:
        if Y < 0.0:
            return np.arctan2(Y, X)+2*np.pi
        else:
            return np.arctan2(Y, X)

def sph_to_rect(r, theta, phi):
    #theta is polar
    #phi is azimuthal
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return np.array([x, y, z])

def rect_to_sph(x, y, z):
    rho = np.sqrt(x**2 + y**2 + z**2)
    phi = theta_cylindrical(x, y) # the azimuthal angle
    theta = np.arccos(z/rho) # the polar angle
    return rho, theta, phi

# return date/time, in format 150102_114501 if date is jan 02 2015, at 11:45:01 am.
def get_time_std():
    millis = str(datetime.datetime.now())[20:22]
    return time.strftime("%Y%m%d-%H%M%S") + millis

def speed_from_KE(KE, m):
    rest = m*c**2
    return c*np.sqrt(1 - (rest/(KE+rest))**2)