from __future__ import division, absolute_import, print_function, unicode_literals
import time
import numpy as np
import sys
import random
import datetime

# Tracking Magnetic Monopoles Through the Galaxy:
# input starting conditions into the program and it will track the magnetic
# monopole's path through the Bfield given by "A New Model of the Galatic 
# Magnetic Field" Jansson, Farrar (2012).

# units used are all SI units except: distance is in kpc, and angles are in
# degrees, magnetic field strength starts in microGauss but output
# of field function converts it to Tesla in get_B() for use in acceleration.

# variable names generally match those in the Bfield paper.

# INPUT PARAMETERS:
    # E0      = sys.argv[1]
    # theta0  = sys.argv[2]     # the initial direction the MP is moving
    # phi0    = sys.argv[3]
    # x0      = sys.argv[4]
    # y0      = sys.argv[5]
    # z0      = sys.argv[6]

# META PARAMETERS
writeToFile = False     # determines if run data is written to a file.
                        # it is overwritten to "True" for batch runs.
outputfilename = "monopoles"    # default name of the optional output file. 
                                # path and ext will be appended later.

simulatebool = False    # determines if trajectory simulation should occur

diskplotbool   = False  # determines if the various plots are made.
trajectorybool = False
vecplotbool    = True
xycrossbool    = False
xzcrossbool    = False

minRhoBool = True   # if true, record distance from galactic center 
                    # at each step

bfieldtests = False # if true, find avg & min/max of bfield

timedirection = 1   # track forwards in time if 1, backwards if -1

relativistic = True # decides whether to use relativistic or nonrelativistic
                    # expression for acceleration.
                    # nonrelativistic computes much more quickly.

halobool = True     # turn on or off various components of the field
xfieldbool = True   # the random field is not used as it is unimportant 
diskbool = True     # this application

print("Note: Run me with python -i monopoles.py to keep terminal " \
    "open afterwards.")
# =======DEFAULT STARTING CONDITIONS=======
# the dirac magnetic charge of the particle (Ampere*kpc)
q_b = 1.06655979e-28

# mass (kg, enter in GeV). An estimate from Wick (2002).
m = 1e15*1.78266184e-27

# stop the program once rho > stopdistance.
# redo the program  if  rho > stopdistance + stopdistancesensitivity
# with a smaller timestep.
STOPDISTANCE = 20.01
STOPDISTANCE_SENSITIVITY = 0.1

# after traversing an arclength of distance_to_halt (kpc), stop the 
# program even if it hasn't gone past stopdistance.
# Currently disabled.
distance_to_halt = 40

c = 9.715611713408621e-12           # light speed in kpc/s

# magnitude of a vector
def mag(V):
    return np.sqrt(V[0]**2+V[1]**2+V[2]**2)

def sphToRect(r, theta, phi):
    #theta is polar
    #phi is azimuthal
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return np.array([x, y, z])

# starting position (kpc), position of earth is x=-8.33937979kpc
# and default timestep (seconds)
r_earth = np.array([-8.5, 0, 0])
# pos0 = np.array([-19.8, 1.5, 1])
pos0 = r_earth
# v_toEarth = np.array([-8.5,0,0]) - pos0       # point velocity to the earth
# vel0 = c*(.4e-3)*(v_toEarth)*(1.0/mag(v_toEarth))

# E0, enter in GeV.
E0 = (10**11)*1.6827e-49
vel0mag = E0*c/np.sqrt(c**4*m**2+E0**2)
theta0 = 1.7
phi0 = 3.14
vel0 = vel0mag*sphToRect(1, theta0, phi0)
dt_default = 5.0e10

# ==set directory to location of this file==
# import os
# abspath = os.path.abspath("__file__") #path of this file
# print("abspath:",abspath)
# dname = os.path.dirname(abspath)
# print("dname:",dname)
# os.chdir(dname)

# =======CONSTANT PARAMETERS=======
c = 9.715611713408621e-12           # light speed in kpc/s

microGtoT = 1e-10
GeVtokg = 1.78266184e-27
mtokpc = 3.24077929e-20

# radius in the Galaxy after which the B field is negligible (kpc)
r_G = 20
# these are the constant portions of the expressions for the maximum
# stepsizes in x and v.
# evaluate them now so that we don't have to re-evaluate them every time.
xMAXconst = 5/(8*r_G)
vMAXconst = c/(200*r_G)


# =======INPUTTED STARTING CONDITIONS=======
try: 
    # first, the parameters needed for earth-start runs
    E0 = float(sys.argv[1])
    vel0mag = E0*c/np.sqrt(c**4*m**2+E0**2)
    theta0=float(sys.argv[2])
    phi0=float(sys.argv[3])
    vel0 = sphToRect(vel0mag, theta0, phi0)

    # assume it is a batch run now. Change parameters accordingly:
    writeToFile = True      # if true, write run data to a file
    simulatebool = True     # if true, simulate mpole trajectory
    diskplotbool   = False  # if true, this type of plot is made
    trajectorybool = False
    vecplotbool    = False
    xycrossbool    = False
    xzcrossbool    = False
    minRhoBool = True   # if true, record distance from galactic center 
                        # at each step.
    bfieldtests = False # if true, find avg, min & max of bfield
    
    # now, the parameters which only vary in 25kpc-start runs
    pos0 = np.array([float(sys.argv[4]), float(sys.argv[5]), 
        float(sys.argv[6])])
except IndexError:
    print("No / not enough parameters inputted. " \
        "Using default parameter values.")



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


def subregioncheck(r, theta, i1, i2):
    if r > r_i(theta, i1) and r < r_i(theta, i2):
        return True
    else:
        return False

# the following function answers: am I in between regions i1 and i2?
def regioncheck(r, theta, i1, i2):
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

def regiondebug(r, theta): # outputs region number, 0 to 7.
    regionsfound = 0
    region_numbers = False
    
    for n in range(8):
        if regioncheck(r, theta, n, (n+1) % 8):
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

# def Theta_X(R,Z):
#     if abs(Z)>=np.tan(Theta0_X)*np.sqrt(R)-np.tan(Theta0_X)*rc_X:
#         return np.arctan(abs(Z)/(R-r_p(R,Z)))  #think about atan and maybe problems with quadrants
#     else:
#         return Theta0_X

# striation parameters (the striated field is currently unused)
#gamma= 2.92
#alpha= 2.65 #taken from page 9, top right paragraph
#beta= gamma/alpha-1

# =======OTHER PRELIMINARY DEFINITIONS=======

# runge-kutta scalars (if k is a slope, not a step):

# with f = dx/dt and g = dv/dt:

# k1x = f(t,          x,               y                 )
# k1y = g(t,          x,               y                 )

# k2x = f(t + a2*dt,  x + (b21*k1x)*dt,  y + (b21*k1y)*dt)
# k2y = g(t + a2*dt,  x + (b21*k1x)*dt,  y + (b21*k1y)*dt)

# ...

# k6x = f(t + a2*dt,  x + (b61*k1x)*dt + (b62*k2x)*dt + (b63*k3x)*dt+... +b65*k5x*dt, y+..)

# AKA since we have unchanging t,
# k['x'][6] = f(x0 + sum([bij[6][n]*k['x'][n]*dt for n in range(1, 6-1)]), vel0 + sum([bij[6][n]*k['v'][n]*dt for n in range(1, 6-1)]) )    
# k['v'][6] = g(x0 + sum([bij[6][n]*k['x'][n]*dt for n in range(1, 6-1)]), vel0 + sum([bij[6][n]*k['v'][n]*dt for n in range(1, 6-1)]) )

# ...but g is accfunction(x,v) and f is just v (which is just a function of v, not x, v). so we have

# k['x'][6] = vel0 + sum([bij[6][n]*k['v'][n]*dt for n in range(1, 6 - 1)])
# k['v'][6] = accfunction(x0 + sum([bij[6][n]*k['x'][n]*dt for n in range(1, 6 - 1)]), vel0 + sum([bij[6][n]*k['v'][n]*dt for n in range(1, 6 - 1)]) )

# f(x0 + sum([bij[6][n]*k['x'][n]))


# Parameters for adaptive timestep Runge-Kutta 5th order.
# ai is unused because of time independence of the acceleration & velocity ODEs.
# ai = np.array([0.0,        0.0, 0.2,         0.3,         0.6,       1.0,    0.875])
# bij = np.array([np.zeros(6),
#                np.array([0.,        0.2,      0.,        0.,           0.,       0.]),
#                np.array([0.,      0.075,   0.225,        0.,           0.,       0.]),
#                np.array([0.,        0.3,    -0.9,       1.2,           0.,       0.]),
#                np.array([0.,     -11/54,     2.5,    -70/27,        35/27,       0.]),
#                np.array([0., 1631/55296, 175/512, 575/13824, 44275/110592, 253/4096]) ])
               
# ci     = np.array([0.0,     37/378, 0.0,     250/621,     125/594,       0.0, 512/1771])
# cistar = np.array([0.0, 2825/27648, 0.0, 18575/48384, 13525/55296, 277/14336,     0.25])


# return date/time, in format 150102_114501 if date is jan 02 2015, at 11:45:01 am.
def gettimestd():
    millis = str(datetime.datetime.now())[20:22]
    return time.strftime("%Y%m%d-%H%M%S")+millis

def exitstatus(pos):
    if mag(pos) >= STOPDISTANCE:
        return True

# rectangular coordinates to polar
def Theta(X, Y): # this theta goes from 0 to 2pi
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

# The lorentz factor
def Gam(vel):
    magvel = mag(vel)
    return c/(np.sqrt(c + magvel)*np.sqrt(c-magvel))

# Kinetic Energy
def KE(vel):
    if Gam(vel)-1.0000001 > 0:
        return m*c**2*(Gam(vel)-1)*5.942795e48     # KE, converted to GeV
    else:
        return m*mag(vel)**2*5.942795e48/2

# xfield function as defined in Farrar
def halofield(pos, r, phi_hat):
    if pos[2] >= 0.0:
        return np.exp(-abs(pos[2])/z_0)*L(pos[2], 
            h_disk, w_disk)* B_n*(1-L(r, r_n, w_h)) * phi_hat

    else:
        return np.exp(-abs(pos[2])/z_0)*L(pos[2], h_disk, w_disk)* \
        B_s*(1-L(r, r_s, w_h)) * phi_hat

def xfield(pos, r, r_p, b_X ):          
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

def diskfield(pos, r, theta, phi_hat):
    if r >= 3.0 and r < 5.0:
        return b_ring*(1.0-L(pos[2], h_disk, w_disk))*phi_hat
            
    if r >= 5.0:
        diskhat = np.array([np.sin(i-theta), np.cos(i-theta), 0.0])
        if regioncheck(r, theta, 7, 0):    # region 1
            return (b1/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 0, 1):  # region 2
            return (b2/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 1, 2):  # region 3
            return (b3/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 2, 3):  # region 4
            return (b4/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 3, 4):  # region 5
            return (b5/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 4, 5):  # region 6
            return (b6/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 5, 6):  # region 7
            return (b7/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
        elif regioncheck(r, theta, 6, 7):  # region 8
            return (b8/r)*(1.0-L(pos[2], h_disk, w_disk))*diskhat
        
    else:
        return [0.0, 0.0, 0.0]

def get_B(pos):
    r = np.sqrt(pos[0]**2.0+pos[1]**2.0)
    bfield = np.array([0.0, 0.0, 0.0])  # B is 0 for r>1, r<20
    
    if mag(pos) > 1.0 and r < 20.0 and r != 0:
        theta = Theta(pos[0], pos[1])
        r_p = R_p(r, pos[2])
        b_X = b_Xf(r_p)
        phi_hat = np.array([-pos[1]/r, pos[0]/r, 0.0])

        #can be optimized by removing if statements
        if halobool:
            bfield += halofield(pos, r, phi_hat)
        if xfieldbool:
            bfield += xfield(pos, r, r_p, b_X)
        if diskbool:
            bfield += diskfield(pos, r, theta, phi_hat)
        
        bfield = bfield*microGtoT  #convert to Tesla
    return bfield

# FIND THE MAXIMUM B MAGNITUDE ON An XY CROSS SECTION (debug)
#maxB=0.0
#for xxx in np.linspace(-20,20,num=200):
#    for yyy in np.linspace(-20,20,num=200):
#        if xxx**2+yyy**2<25:
#            if mag(get_B([xxx,yyy,.01]))>maxB:
#                maxB=mag(get_B([xxx,yyy,.01]))
#print("MAXB IS",maxB)

#returns relativistic acceleration. only necessary for large velocity.
def acc_relativistic(pos, vel):   
    bfield = get_B(pos)
    gam = Gam(vel)
    force = q_b*bfield    # In kg*kpc/s^2
    acc = (1.0/(gam*m*c**2)) * (c**2*force - np.dot(force, vel)*vel)
    return acc

def vel(pos, vel):
    return vel

#this is really only a function of position.
def acc_classical(pos, vel):
    bfield = get_B(pos)
    acc  = bfield*(q_b/m)
    return acc

if relativistic: accfunction = acc_relativistic
else:            accfunction = acc_classical

#k[i][j]  is the ith k of the j=0 (position) or j=1 (velocity), it is the slope
#slope=[0, [0,0], [0,0], [0,0], [0,0], [0,0], [0,0] ]
#k = [0, {'x':0,'v':0}, {'x':0,'v':0}, {'x':0,'v':0}, {'x':0,'v':0}, {'x':0,'v':0}, {'x':0,'v':0}]
k = {'x' : ['null', 'k1', 'k2', 'k3', 'k4', 'k5', 'k6'], 'v' : ['null', 'k1', 'k2', 'k3', 'k4', 'k5', 'k6']}

#okay, rk5. f1(t,x) gives the velocity at x, t. 
#we have 2 relationships: 1. x''=acc(pos,vel) 2. x'=acc*t
#           f2(t,x) gives the acceleration at x,t.
#           f2 is a plug in from accfunction.
#           f1 is dx/dt=
def eulerStep(pos, vel, drdt, dvdt, dt):
    dr = drdt(pos,vel)*dt
    dv = dvdt(pos,vel)*dt
    return dr, dv
    
def rk4step(pos, vel, drdt, dvdt, dt):
    #return the step to be added, NOT the updated vector, so that dr can be 
    #added to distance_traveled.
    k1 = np.array([vel*dt,                    accfunction(pos,           vel)*dt           ])
    k2 = np.array([(vel+k1[1]/2.0)*dt,        accfunction(pos+k1[0]/2.0, vel+k1[1]/2.0)*dt ])
    k3 = np.array([(vel+k2[1]/2.0)*dt,        accfunction(pos+k2[0]/2.0, vel+k2[1]/2.0)*dt ])
    k4 = np.array([(vel+k3[1])*dt,            accfunction(pos+k3[0],     vel+k3[1])*dt     ])
    
    pos_step = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
    # tested: mutability of np array does not mess this up (changing vel, pos does not chg k1)
    vel_step = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0
    
    return pos_step, vel_step
    
    
def rk5Step(xn, vn, dt):
    #adaptive timestep 5th order method (in progress)
    # set the first k's to build up the others with
    k['x'][1] = vn                              # velocity k1
    k['v'][1] = accfunction(xn, vn)         # position k1
    
    # set the rest of the k's (see numerical recipes in c pg 714)
    for i in range(2, 7):
        # summing a series of numpy arrays with sum() DOES result in their vector sum.
        # CHANGE THIS TO NUMPY SUM IF SUCH A THING EXISTS, FOR OPTIMIZATION
        k['x'][i] = vn + sum([bij[i][j]*k['v'][j] for j in range(1, i)])*dt
        k['v'][i] = accfunction(xn + sum([bij[i][j]*k['x'][j] for j in range(1, i)])*dt, vn + sum([bij[i][j]*k['v'][j] for j in range(1, i)])*dt)
    
#   5th order method:
    xn_O5 = xn + sum([ci[i]*k['x'][i]*dt for i in range(1, 7)])
    vn_O5 = xn + sum([ci[i]*k['v'][i]*dt for i in range(1, 7)])
    
#    4th order method:
#    xn_O4 = xn + sum([cistar[i]*k['x'][i]*dt for i in range(1, 7)])
#    vn_O4 = xn + sum([cistar[i]*k['v'][i]*dt for i in range(1, 7)])
    
#   error.f =5thorder-(4thorder)
    Deltax = xn_O5 - (xn + sum([cistar[i]*k['x'][i]*dt for i in range(1, 7)]) )
    Deltav = vn_O5 - (xn + sum([cistar[i]*k['v'][i]*dt for i in range(1, 7)]) )
    
    # our maximum errors before we redo a step.
    Deltax0 = xMAXconst*distance_tracked/iterations
    Deltavel0 = vMAXconst*distance_tracked/iterations
    
    # the step size we SHOULD have taken to get the right error
    dtgood = dt*np.pow(abs(DeltaxMAX), 0.2)
    
    # if the largest acceptable error is bigger than our approximated error:
    if (DeltaxMAX > max(Deltax)) and (DeltavMAX > max(Deltav)): # THIS IS BAD BC DELTAX IS A VECTOR
        # our step was good. return the updated x and v, and a new stepsize dt for next time.
        xn_O5 = xn + sum([ci[i]*k['x'][i]*dt for i in range(1, 7)])
        vn_O5 = xn + sum([ci[i]*k['v'][i]*dt for i in range(1, 7)])
        
        dt_next = dt*(Drjolrkefkelqfaoshjolasjd)
        return [xn_O5, vn_O5]
    else:
        print("This part isn't done yet!!")
        
    
    
    
    
    
# idea: make a relativistic and nonrel paradigm. when the speed is nonrel (difference
# would be less than machine precision) relativistic=False. Every X iterations, check if 
# the paradigm should be changed, where X is the least number of iterations it would take
# to change the velocity by, say, 1% of the speed of light.


# =======PROGRAM BEGINS=======
def runTrajectory(pos0,vel0, dt):
    # ======== PART 1: PREPARATION ========    
    #==optional parameters & other stuff==
    tstartsecs = time.time()
    maxvelocity = mag(vel0)
    iterations = 0
    arclength_in_rho1kpc_region = 0.0
    #for writing to a file
    if (writeToFile == True):
        tstartstd = gettimestd()
        a = "data/"+outputfilename+gettimestd()+".txt"
        f = open(a, 'a')
        f.write("=====STARTING CONDITIONS=====\n") 
        f.write("start time (YYYYMMDD-hhmmssms)\n")
        f.write(gettimestd()+"\n")
        f.write("vel0 (kpc/s):\n")
        f.write(str(vel0)+"\n")
        f.write("mag(vel0):\n")
        f.write(str(mag(vel0))+"\n")
        f.write("mag(vel0)/c\n")
        f.write(str(mag(vel0)/c)+"\n")
        f.write("E0 (kg*kpc^2/s^2):\n")
        f.write(str(E0)+"\n")
        f.write("theta0 (rads):\n")
        f.write(str(theta0)+"\n")
        f.write("phi0 (rads):\n")
        f.write(str(phi0)+"\n")
        f.write("q_b\n")
        f.write(str(q_b)+"\n")
        f.write("m (kg)\n")
        f.write(str(m)+"\n")
        f.write("dt (s)\n")
        f.write(str(dt)+"\n")
        f.write("distance to halt (kpc)\n")
        f.write(str(distance_to_halt)+"\n")     # line 22

    if exitstatus(pos0) == True:
        print("ERROR: Started already outside STOPDISTANCE.")
        print("Adding to list of failed runs.")
        f.write("ERROR: RUN FAILED. ERROR CODE:\n")     # line 23
        f.write("POS0ERROR")                            # line 24
        sys.exit(1)

    #for graphing
    xlist, ylist, zlist = [np.array([]) for _ in range(3)]
    rholist = np.array([])

    # ======== PART 2: LOOP ========
    pos = np.copy(pos0)
    vel = np.copy(vel0)
    acc = accfunction(pos, vel) # python complains unless acc is defined 
                                # outside the loop.
    distance_tracked = 0.0   #halt the loop if the monopole travels too far.
    
    # if v>c then return vel=(c,c,c), which is error code for vel > c
    if mag(vel) >= c: 
        print("mag(vel)>c. Returning (c,c,c) error.")
        return pos, np.array([c, c, c])

    print("Starting traversal. Time is moving ",end="")
    if timedirection==1:  print("FORWARDS.")
    if timedirection==-1: print("BACKWARDS.")

    if relativistic == True:  print("Using RELATIVISTIC acceleration.")
    if relativistic == False: print("Using CLASSICAL acceleration.")

    while mag(pos) < STOPDISTANCE:
        
        iterations += 1
        
        # update trajectory. Use timedirection =-1 to track backwards.
        acc = accfunction(pos, vel)
        pos_step, vel_step = rk4step(pos, vel, vel, acc, dt)
        pos += timedirection*pos_step
        vel += timedirection*vel_step
        distance_tracked += mag(pos_step)

        # if v>c then return vel=(c,c,c), which is error code for vel > c
        if mag(vel) >= c: return pos, np.array([c, c, c])
        
        # FOR GRAPHING:
        if trajectorybool:
            xlist = np.concatenate((xlist, np.array([pos[0]])))
            ylist = np.concatenate((ylist, np.array([pos[1]])))
            zlist = np.concatenate((zlist, np.array([pos[2]])))

        if minRhoBool:
            rholist = np.concatenate((rholist, np.array([mag(pos)])))        
        
        if iterations % 10000 == 0:
            bf = get_B(pos)
            disp = pos - pos0
            print("============done with iteration no.", iterations, "============")
            print("pos-pos0 = [%.3g, %.3g, %.3g]. |pos-pos0| = %.3gkpc" % (disp[0], disp[1], disp[2],mag(disp)))
            print("vel/c = [%.3g, %.3g, %.3g]. |vel|/c = %.3g" % (vel[0]/c, vel[1]/c, vel[2]/c,mag(vel)/c))            
            print("acc/c = [%.3g, %.3g, %.3g]. |acc|/c = %.3g/s" % (acc[0]/c, acc[1]/c, acc[2]/c,mag(acc)/c))            
            print("bfield = [%.3g, %.3g, %.3g]. |bfield| = %.3gT" % (bf[0], bf[1], bf[2],mag(bf)))
            print("arc length = %.3gkpc" % distance_tracked)
            print("simulated time = %.3gs" % (iterations*timedirection*dt))
            print("real runtime = %.3g real seconds" % (time.time()-tstartsecs))
            print()
            if relativistic == False and mag(vel) > 0.2*c:
                print("\n=====================================================")
                print("====================ERROR===========================")
                print("=== mag(velocity) > 0.2c but you are not using    ===")
                print("===      the relativistic acceleration!!!         ===")
                sys.exit(1)
            
        if mag(pos) < 1:
            arclength_in_rho1kpc_region += mag(pos_step)
    
    # print some stuff after the run?
    # print("iterations:", iterations)
    # print("Final pos:",pos)
    # print("Final vel/c:\n",vel/c)
    # print()
    
    if exitstatus(pos) == True:
        print("Stopped because mag(pos) = %.3f >= %.2f." % 
            (mag(pos),STOPDISTANCE))
    else:
        print("Stopped because distance traveled >",distance_to_halt)
    if minRhoBool:
        print("minimum rho=mag(pos) was %.3g." % min(rholist))

    distance_from_start = mag(pos-pos0)
    theta_f = np.arccos(pos[2]/mag(pos))
    phi_f = Theta(pos[0], pos[1])
    
    if (writeToFile == True):
        f.write("=====FINAL CONDITIONS=====\n")
        f.write("finish time (YYYYMMDD-hhmmssms)\n")
        f.write(gettimestd()+"\n")
        f.write("vel_final (kpc/s):\n")
        f.write(str(vel)+"\n")
        f.write("mag(vel_final) (kpc/s)\n")
        f.write(str(mag(vel))+"\n")
        f.write("mag(vel_final)/c\n")
        f.write(str(mag(vel)/c)+"\n")
        f.write("pos_final (kpc)\n")
        f.write(str(pos)+"\n")
        f.write("mag(pos_final) (kpc)\n")
        f.write(str(mag(pos))+"\n")
        f.write("theta_final (rads)\n")
        f.write(str(theta_f)+"\n")
        f.write("phi_final (rads)\n")
        f.write(str(phi_f)+"\n")
        f.write("distance tracked (kpc)\n")
        f.write(str(distance_tracked)+"\n")
        f.write("distance from start (kpc)\n")
        f.write(str(distance_from_start)+"\n")
        f.write("Kinetic Energy (kg*kpc^2/s^2)\n")
        f.write(str(KE(vel))+"\n")
        f.write("maxvelocity (kpc/s)\n")
        f.write(str(maxvelocity)+"\n")
        f.write("maxvelocity/c\n")
        f.write(str(maxvelocity/c)+"\n")
        f.write("arclength_in_rho1kpc_region (kpc)\n")
        f.write(str(arclength_in_rho1kpc_region)+"\n")
        f.write("time after t=0 start (s):\n")
        f.write(str(iterations*dt*timedirection)+"\n")
        f.write("iterations\n")
        f.write(str(iterations)+"\n")
        f.write("real runtime (real s)\n")
        f.write(str(time.time()-tstartsecs)+"\n")
        f.write("acc_final (kpc/s/s)\n")
        f.write(str(acc)+"\n\n\n")
        f.write("exit status\n")
        f.write(str(exitstatus(pos)))
        f.close()

    if trajectorybool:
        with open("grapher.py") as f:
            code = compile(f.read(), "grapher.py", 'exec')
            print()
            exec(code)
    return pos,vel

if not halobool:
    print("WARNING: The halo field component is disabled (halobool==False)")
if not xfieldbool:
    print("WARNING: The x field component is disabled (xfieldbool==False)")
if not diskbool:
    print("WARNING: The disk field component is disabled (diskbool==False)")
if not simulatebool:
    if trajectorybool:
        print("ERROR: trajectorybool==True but simulatebool==False ! Halting.\n")
        sys.exit(0)
    if minRhoBool:
        print("")
        print("ERROR: minRhoBool==True but simulatebool==False ! Making minRhoBool==False.\n")
        minRhoBool=False
    print("WARNING: You have chosen not to simulate the trajectory. (simulatebool==False)")

if simulatebool:
    pos_final, vel_final = runTrajectory(pos0, vel0, dt_default)

    runcount_finalposerror = 0
    while( mag(pos_final) > STOPDISTANCE+STOPDISTANCE_SENSITIVITY):
        print("ERROR: Must redo run for final position accuracy.")
        print("mag(pos_final) was %.3fkpc, which is greater than the " \
            "threshhold of %.3f + %.3f = %.3fkpc." %
            (mag(pos_final), STOPDISTANCE, STOPDISTANCE_SENSITIVITY,
                STOPDISTANCE+STOPDISTANCE_SENSITIVITY))
        print("Starting redo attempt number %d, with a halved timestep.\n" % 
            (runcount_finalposerror+1))
        dt_default = dt_default/2
        pos_final, vel_final = runTrajectory(pos0, vel0, dt_default)
        runcount_finalposerror+=1

    runcount_speederror = 0
    while( mag(vel_final)/mag(np.array([c, c, c])) > 0.9):
        # returning vel=c, c, c is error code for mag(vel)>c. 
        # redo run with halved timestep
        print("====HAD TO REDO RUN BECAUSE mag(vel) > c ====")
        print("Starting redo attempt number %d with a halved timestep." % 
            (runcount_speederror+1))
        dt_default = dt_default/2
        pos_final, vel_final = runTrajectory(pos0, vel0, dt_default)
        runcount_speederror += 1

if bfieldtests:
    with open("bfieldtests.py") as f:
        code = compile(f.read(), "bfieldtests.py", "exec")
        print()
        exec(code)

if (diskplotbool or vecplotbool or xycrossbool or xzcrossbool) and not trajectorybool:
    # trajectory was already done in the runTrajectory function,
    # as it needed local variables. All other plots were also done at the time.
    # Only do this if they haven't been done yet.
    with open("grapher.py") as f:
        code = compile(f.read(), "grapher.py", 'exec')
        print()
        exec(code)
print()