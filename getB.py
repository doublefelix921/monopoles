from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass
try: range = xrange
except NameError: pass

import numpy as np
import usefulfunctions as f

# META PARAMETERS: LEAVE DEFAULT TRUE!!!
HaloFieldBool = True     # turn on or off various components of the field
XFieldBool = True   # the random field is not used as it is unimportant 
DiskFieldBool = True     # this application

microG_TO_T = f.microG_TO_T

# magnitude of a vector
mag = f.mag

# rectangular coordinates to polar
theta_cylindrical= f.theta_cylindrical

# =======DISK PARAMETERS, WITH B IN microGauss. =======
# =======(Converted to Tesla when field is calculated)=======
b1 = 0.1
b2 = 3.0
b3 = -0.9
b4 = -0.8
b5 = -2.0
b6 = -4.2
b7 = 0.0
b8 = 2.7
b_ring = 0.1
h_disk = 0.4
w_disk = 0.27

# =======HALO PARAMETERS & FUNCTIONS=======
B_n = 1.4
B_s = -1.1
r_n = 9.22
r_s = 16.7  # NOTE: parameter has a large error.
w_h = 0.2
z_0 = 5.3
# this is the opening "pitch" angle of the logarithmic spiral boundaries.
i = np.radians(11.5)
# these are the radii at which the spirals cross the x-axis
r_negx = np.array([5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5])

# halo spiral boundary function
def r_i(T, I):
    return r_negx[I]*np.exp((T-np.pi)*np.tan(i))


def subregion_check(r, theta, i1, i2):
    if r > r_i(theta, i1) and r < r_i(theta, i2):
        return True
    else:
        return False

# the following function answers: am I in between regions i1 and i2?
def region_check(r, theta, i1, i2):
    if i1 != 7:
        if r > r_i(theta-2*np.pi, i1) and r < r_i(theta-2*np.pi, i2):
            return True
        if r > r_i(theta, i1) and r < r_i(theta, i2):
            return True
        if r > r_i(theta+2*np.pi, i1) and r < r_i(theta+2*np.pi, i2):
            return True
        else:
            return False
    if i1 == 7:
        if r > r_i(theta-2*np.pi, i1) and r < r_i(theta, i2):
            return True
        if r > r_i(theta, i1) and r < r_i(theta+2*np.pi, i2):
            return True
        if r > r_i(theta+2*np.pi, i1) and r < r_i(theta+4*np.pi, i2):
            return True
        else:
            return False

def L(Z, H, W):          # disk-halo transition function
    return 1/(1.0+np.exp(-2.0*(abs(Z)-H)/W))

def region_debug(r, theta): # outputs region number, 0 to 7.
    regionsfound = 0
    region_numbers = False
    
    for n in range(8):
        if region_check(r, theta, n, (n+1) % 8):
            regionsfound += 1
            region_numbers = (n, (n+1) % 8)
    return regionsfound, region_numbers

# =======X-FIELD PARAMATERS & FUNCTIONS=======
B_X = 4.6
Theta0_X = np.radians(49.0)
rc_X = 4.8
r_X = 2.9

def R_p(R, Z):
    if abs(Z) >= np.tan(Theta0_X)*(R-rc_X):
        return R*rc_X/(rc_X+abs(Z)/np.tan(Theta0_X))

    else:
        return R-abs(Z)/np.tan(Theta0_X)


def b_Xf(R_P):
    return B_X*np.exp(-R_P/r_X)

# xfield function as defined in Farrar
def halo_field(pos, r, phi_hat):
    if pos[2] >= 0.0:
        return np.exp(-abs(pos[2])/z_0)*L(pos[2], 
            h_disk, w_disk)* B_n*(1-L(r, r_n, w_h)) * phi_hat

    else:
        return np.exp(-abs(pos[2])/z_0)*L(pos[2], h_disk, w_disk)* \
        B_s*(1-L(r, r_s, w_h)) * phi_hat

def x_field(pos, r, r_p, b_X ):
    if pos[2] != 0:
        bhat_X = 1/np.sqrt(pos[2]**2+(r-r_p)**2)*np.array([pos[0]*(r-r_p)/r, 
            pos[1]*(r-r_p)/r, pos[2]])
        if pos[2] < 0:
            bhat_X *= -1
        if abs(r_p) < rc_X:
            return b_X*(r_p/r)**2*bhat_X

        if abs(r_p) >= rc_X:
            return b_X*(r_p/r)*bhat_X
    else:
        return b_X*np.array([0, 0, 1])

def disk_field(pos, r, theta, phi_hat):
    if r >= 3.0 and r < 5.0:
        return b_ring*(1.0-L(pos[2], h_disk, w_disk))*phi_hat
            
    if r >= 5.0:
        diskhat = np.array([np.sin(i-theta), np.cos(i-theta), 0.0])
        if region_check(r, theta, 7, 0):    # region 1
            return (b1/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 0, 1):  # region 2
            return (b2/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 1, 2):  # region 3
            return (b3/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 2, 3):  # region 4
            return (b4/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 3, 4):  # region 5
            return (b5/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 4, 5):  # region 6
            return (b6/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 5, 6):  # region 7
            return (b7/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif region_check(r, theta, 6, 7):  # region 8
            return (b8/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
    else:
        return [0.0, 0.0, 0.0]

def get_B(pos):
    r = np.sqrt(pos[0]**2.0+pos[1]**2.0)
    bfield = np.array([0.0, 0.0, 0.0])  # B is 0 for r>1, r<20
    
    if mag(pos) > 1.0 and r < 20.0 and r != 0:
        theta = theta_cylindrical(pos[0], pos[1])
        r_p = R_p(r, pos[2])
        b_X = b_Xf(r_p)
        phi_hat = np.array([-pos[1]/r, pos[0]/r, 0.0])

        #can be optimized by removing if statements
        if HaloFieldBool:
            bfield += halo_field(pos, r, phi_hat)
        if XFieldBool:
            bfield += x_field(pos, r, r_p, b_X)
        if DiskFieldBool:
            bfield += disk_field(pos, r, theta, phi_hat)
        
        bfield = bfield*microG_TO_T  #convert to Tesla
    return bfield