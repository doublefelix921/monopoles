from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass

import traceback
import time
import numpy as np
import sys
import random
import datetime
import code # use code.interact(local=locals()) 
            # to return to command line at that point for debugging
import getB as gb

# Tracking Magnetic Monopoles Through the Galaxy:
# input starting conditions into the program and it will track the magnetic
# monopole's path through the Bfield given by "A New Model of the Galatic 
# Magnetic Field" Jansson, Farrar (2012).

# units used are all SI units except: distance is in kpc, and angles are in
# degrees, magnetic field strength inputs are in microGauss but calculations
# are done in Tesla so the gb.get_B() function does conversions before output.

# variable names generally match those in the Bfield paper.

# INPUT PARAMETERS:
    # E0      = sys.argv[1]
    # vtheta0  = sys.argv[2]     # the initial direction the MP is moving
    # vphi0    = sys.argv[3]
    # x0      = sys.argv[4]
    # y0      = sys.argv[5]
    # z0      = sys.argv[6]

# META PARAMETERS
writeToFile = True     # determines if run data is written to a file.
                        # it is overwritten to "True" for batch runs.
outputfilename = "monopoles"    # default name of the optional output file. 
                                # path and ext will be appended later.

simulatebool = True    # determines if trajectory simulation should occur

diskplotbool   = False  # determines if the various plots are made.
trajectorybool = False
vecplotbool    = False
xycrossbool    = False
xzcrossbool    = False

minRhoBool = True   # if true, record distance from galactic center 
                    # at each step

bfieldtests = False # if true, find avg & min/max of bfield

timedirection = 1   # track forwards in time if 1, backwards if -1

relativistic = True # decides whether to use relativistic or nonrelativistic
                    # expression for acceleration.
                    # nonrelativistic computes much more quickly.

gb.halobool = True     # turn on or off various components of the field
gb.xfieldbool = True   # the random field is not used as it is unimportant 
gb.diskbool = True     # this application

# =======DEFAULT STARTING CONDITIONS=======
# the dirac magnetic charge of the particle (Ampere*kpc?)
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
# If distance_to_halt is -1, don't stop based on distance req
distance_to_halt = -1

c = 9.715611713408621e-12           # light speed in kpc/s
kgkpc2s2_To_GeV = 5.9427943e48

# magnitude of a vector
def mag(V):
    return np.sqrt(V[0]**2+V[1]**2+V[2]**2)

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

def sphToRect(r, theta, phi):
    #theta is polar
    #phi is azimuthal
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return np.array([x, y, z])

def rectToSph(x, y, z):
    rho = np.sqrt(x**2 + y**2 + z**2)
    phi = Theta(x, y) # the azimuthal angle
    theta = np.arccos(z/rho) # the polar angle
    return rho, theta, phi


# starting position (kpc), position of earth is x=-8.33937979kpc
# and default timestep (seconds)
pos_earth = np.array([-8.5, 0, 0])
# pos0 = np.array([-19.8, 1.5, 1])
pos0 = pos_earth
# v_toEarth = np.array([-8.5,0,0]) - pos0       # point velocity to the earth
# vel0 = c*(.4e-3)*(v_toEarth)*(1.0/mag(v_toEarth))

# E0, enter in GeV in the parenthesis
E0 = (10**11)*1.6827e-49
vel0mag = E0*c/np.sqrt(c**4*m**2+E0**2)
vtheta0 = 1.7
vphi0 = 3.14
vel0 = vel0mag*sphToRect(1, vtheta0, vphi0)
dt_default = 5.0e10
dt = dt_default

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

# =======INPUTTED STARTING CONDITIONS=======
try: 
    # first, the parameters needed for earth-start runs
    E0 = float(sys.argv[1])
    vel0mag = E0*c/np.sqrt(c**4*m**2+E0**2)
    vtheta0=float(sys.argv[2])
    vphi0=float(sys.argv[3])
    vel0 = sphToRect(vel0mag, vtheta0, vphi0)

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
    print("No / not enough parameters inputted. ")
    ans = input("Use default parameter values for those which " +
        "weren't inputted and continue?\n y/n ")
    if not (ans.lower() in ["y", "yes", "yy"]):
        sys.exit(1)

# return date/time, in format 150102_114501 if date is jan 02 2015, at 11:45:01 am.
def gettimestd():
    millis = str(datetime.datetime.now())[20:22]
    return time.strftime("%Y%m%d-%H%M%S") + millis

def exitstatus(pos):
    if mag(pos) >= STOPDISTANCE:
        return True

# The lorentz factor
def Gam(vel):
    magvel = mag(vel)
    return c/(np.sqrt(c+magvel)*np.sqrt(c-magvel))

# Kinetic Energy
def KE(vel):
    #return in kg*kpc^2/s^2
    if Gam(vel) > 1.000000001:
        return m*c**2*(Gam(vel)-1)
    else:
        return m*mag(vel)**2/2

# FIND THE MAXIMUM B MAGNITUDE ON An XY CROSS SECTION (debug)
#maxB=0.0
#for xxx in np.linspace(-20,20,num=200):
#    for yyy in np.linspace(-20,20,num=200):
#        if xxx**2+yyy**2<25:
#            if mag(gb.get_B([xxx,yyy,.01]))>maxB:
#                maxB=mag(gb.get_B([xxx,yyy,.01]))
#print("MAXB IS",maxB)

#returns relativistic acceleration. only necessary for large velocity.
def acc_relativistic(pos, vel):   
    bfield = gb.get_B(pos)
    gam = Gam(vel)
    force = q_b*bfield    # In kg*kpc/s^2
    acc = (1.0/(gam*m*c**2)) * (c**2*force - np.dot(force, vel)*vel)
    return acc

def vel(pos, vel):
    return vel

#this is really only a function of position.
def acc_classical(pos, vel):
    bfield = gb.get_B(pos)
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

# idea: make a relativistic and nonrel paradigm. when the speed is nonrel (difference
# would be less than machine precision) relativistic=False. Every X iterations, check if 
# the paradigm should be changed, where X is the least number of iterations it would take
# to change the velocity by, say, 1% of the speed of light.


# =======PROGRAM BEGINS=======
def runTrajectory(pos0, vel0, dt, m, q_b):
    # ======== PART 1: PREPARATION ========    
    #==optional parameters & other stuff==
    tstartsecs = time.time()
    maxvelocity = mag(vel0)
    iterations = 0
    arclength_in_rho1kpc_region = 0.0
    
    if exitstatus(pos0) == True:
        raise Exception("POS0ERROR")

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
        if time.time()-tstartsecs > 300:
            #if run has been going on more than 5mins, timeout error.
            raise Exception("SPEEDTIMEOUTERROR")
        print("mag(vel)>c. Returning (c,c,c) error and trying again.")
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
            if time.time()-tstartsecs > 600:
                raise Exception("TIMEOUTERROR")
            bf = gb.get_B(pos)
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

    if exitstatus(pos) == True:
        print("Stopped because mag(pos)=%.3g >= %.3g." % 
            (mag(pos),STOPDISTANCE))
    else:
        print("Stopped because distance traveled >",distance_to_halt)
    if minRhoBool:
        print("minimum rho=mag(pos) was %.3g." % min(rholist))

    # check if change in KE is unusually large:
    # NOT NECESSARY JUST DO IT AFTER RUNS ARE DONE
    # if (KE(vel)-KE(vel0))*kgkpc2s2_To_GeV > 4e12:
    #     ans1 = input("PAUSING: Change in energy was unusually large:\n" +
    #         ("delta KE = %.5gGeV \n" % ((KE(vel)-KE(vel0))*kgkpc2s2_To_GeV)) +
    #         ("path length was %.3gkpc \n " % distance_tracked) +
    #         "starting conditions were:\n" + 
    #         ("pos0: [x,y,z]=[%.8g, %.8g, %.8g]kpc\n" % (pos0[0], pos0[1], pos0[2])) +
    #         ("vel0: [x,y,z]=[%.8g, %.8g, %.8g]kpc/s\n" % (vel0[0], vel0[1], vel0[2])) +
    #         ("timedirection: %d\n" % timedirection) +
    #         "continue? y/n ")
    #     if not(ans1.lower() in ['y', 'yes', 'yy']):
    #         ans2 = input("are you sure you want to quit? y/n ")
    #         if ans2 in ['y', 'yes', 'yy']:
    #             print("got it. exiting.\n")
    #             sys.exit(1)
    #     else:    

    distance_from_start = mag(pos-pos0)
    thetaF = np.arccos(pos[2]/mag(pos))
    phiF = Theta(pos[0], pos[1])

    if trajectorybool:
        with open("grapher.py") as f:
            code = compile(f.read(), "grapher.py", 'exec')
            print()
            exec(code)

    return pos,vel, distance_tracked, maxvelocity, \
        arclength_in_rho1kpc_region, iterations

def runTrajectoryWrapper(pos0, vel0, dt, m, q_b):
    try:

        posF, velF, distance_tracked, maxvelocity, \
            arclength_in_rho1kpc_region, iterations = \
            runTrajectory(pos0,vel0, dt, m, q_b)
    except:
        print("RUN FAILED, MAKING ERROR FILE!!")
        tb = traceback.format_exc()
        print(tb)
        
        with open( (',ERROR'+gettimestd()), 'a' ) as f2:
            f2.write("pos0:\n")
            f2.write(str(pos0)+'\n')
            f2.write("vel0:\n")
            f2.write(str(vel0)+'\n')
            tb = traceback.format_exc()
            f2.write(tb+"\n")
        time.sleep(3) #let the error msg be seen
        sys.exit(1)
    return posF, velF, distance_tracked, maxvelocity, \
            arclength_in_rho1kpc_region, iterations

# ok now calling functions.
print("Note: Run me with python -i monopoles.py to keep terminal " \
    "open afterwards.")

if not gb.halobool:
    print("WARNING: The halo field component is disabled (halobool==False)")
if not gb.xfieldbool:
    print("WARNING: The x field component is disabled (xfieldbool==False)")
if not gb.diskbool:
    print("WARNING: The disk field component is disabled (diskbool==False)")
if not simulatebool:
    if trajectorybool:
        print("ERROR: trajectorybool==True but simulatebool==False ! Halting.\n")
        sys.exit(0)
    if minRhoBool:
        print("")
        print("ERROR: minRhoBool==True but simulatebool==False ! " +
            "Making minRhoBool==False.\n")
        minRhoBool=False
    print("WARNING: You have chosen not to simulate the trajectory. " +
        "(simulatebool==False)")

if simulatebool:
    tstartstd = gettimestd()
    tstartsecs = time.time()
    posF, velF, distance_tracked, maxvelocity, \
        arclength_in_rho1kpc_region, iterations = \
        runTrajectoryWrapper(pos0, vel0, dt, m, q_b)
    
    runcount_finalposerror = 0
    while( mag(posF) > STOPDISTANCE+STOPDISTANCE_SENSITIVITY):
        print("ERROR: Must redo run for final position accuracy.")
        print("mag(posF) was %.3gkpc, which is greater than the " \
            "threshhold of %.3g + %.3g = %.3gkpc." %
            (mag(posF), STOPDISTANCE, STOPDISTANCE_SENSITIVITY,
                STOPDISTANCE+STOPDISTANCE_SENSITIVITY))
        print("Starting redo attempt number %d, with a halved timestep.\n" % 
            (runcount_finalposerror+1))
        dt = dt/2
        tstartstd = gettimestd()
        tstartsecs = time.time()
        posF, velF, distance_tracked, maxvelocity, \
            arclength_in_rho1kpc_region, iterations = \
            runTrajectoryWrapper(pos0, vel0, dt, m, q_b)
        runcount_finalposerror += 1

    runcount_speederror = 0
    while( mag(velF)/mag(np.array([c, c, c])) > 0.9):
        # returning vel=c, c, c is error code for mag(vel)>c. 
        # redo run with halved timestep
        print("====HAD TO REDO RUN BECAUSE mag(vel) > c ====")
        print("Starting redo attempt number %d with a halved timestep." % 
            (runcount_speederror+1))
        dt = dt/2
        tstartstd = gettimestd()
        posF, velF, distance_tracked, maxvelocity, \
            arclength_in_rho1kpc_region, iterations = \
            runTrajectoryWrapper(pos0, vel0, dt, m, q_b)
        runcount_speederror += 1

    # if you are up to here, the run has succeeded.
    reruns = runcount_finalposerror + runcount_speederror

    # print some stuff after the run?
    print("Run complete with no errors after %d reruns." % reruns)
    print("iterations:", iterations)
    print("posF:   ", posF)
    print("velF/c: ", velF/c)
    print("Initial KE(GeV): %g" % (KE(vel0)*kgkpc2s2_To_GeV))
    print("Final KE(GeV): %g" % (KE(velF)*kgkpc2s2_To_GeV))

    # write results to file
    if writeToFile == True:
        rhoF, thetaF, phiF = rectToSph(posF[0], posF[1], posF[2])
        colheaders = ["start time (YYYYMMDD-hhmmssms)", "v0 (kpc/s)", 
    "mag(v0) (kpc/s)", "mag(v0/c)", "E0 (kg*kpc^2/s^2)", "theta0 (rads)", 
    "phi0 (rads)", "q_b (A*kpc)", "mass (kg)", "dt (s)", 
    "distance to halt (kpc)", "end time (YYYYMMDD-hhmmssms)",
    "vel_final (kpc/s)", "mag(vel_final) (kpc/s)", "mag(vel_final/c)", 
    "pos_final(kpc)", "mag(pos_final) (kpc)", "theta_final (rads)", 
    "phi_final (rads)", "distance tracked (kpc)", "distance from start (kpc)",
    "KE_final (kg*kpc^2/s^2)", "max velocity (kpc/s)", "max velocity/c",
    "arclength in rho1kpc region (kpc)", "time after t=0 start (s)",
    "iterations", "real runtime (real s)", "acc_final (kpc/s/s)",
    "exit status (True/False)"]


        outputpath = "data/"+outputfilename+gettimestd()+".txt"
        f = open(outputpath, 'a')

        # initial conditions
        f.write("start time (YYYYMMDD-hhmmssms)\n")
        f.write(tstartstd+"\n")
        f.write("[x0,y0,z0]:\n")
        f.write(str(pos0)+"\n")
        f.write("[vtheta0, vphi0] (kpc,rads,rads):\n")
        f.write(str(np.array([vtheta0, vphi0])) + "\n")
        f.write("vel0 (kpc/s):\n")
        f.write(str(vel0)+"\n")
        f.write("E0 (kg*kpc^2/s^2):\n")
        f.write(str(KE(vel0))+"\n")
        f.write("q_b\n")
        f.write(str(q_b)+"\n")
        f.write("m (kg)\n")
        f.write(str(m)+"\n")
        f.write("dt (s)\n")
        f.write(str(dt)+"\n")
        f.write("distance_to_halt (kpc)\n")
        f.write(str(distance_to_halt) + "\n") 

        # final conditions
        f.write("finish time (YYYYMMDD-hhmmssms)\n")
        f.write(gettimestd()+"\n")
        f.write("[xf, yf, zf] (kpc)\n")
        f.write(str(posF)+"\n")
        f.write("[rhoF, thetaF, phiF] (kpc,rads,rads):\n")
        f.write(str([rhoF, thetaF, phiF]) + "\n")
        f.write("velF (kpc/s):\n")
        f.write(str(velF)+"\n")
        f.write("distance tracked (kpc)\n")
        f.write(str(distance_tracked)+"\n")
        f.write("distance from start (kpc)\n")
        f.write(str(mag(posF-pos0))+"\n")
        f.write("Final KE (kg*kpc^2/s^2)\n")
        f.write(str(KE(velF))+"\n")
        f.write("maxvelocity/c \n")
        f.write(str(maxvelocity/c)+"\n")
        f.write("arclength_in_rho1kpc_region (kpc)\n")
        f.write(str(arclength_in_rho1kpc_region)+"\n")
        f.write("time after t=0 start (s):\n")
        f.write(str(iterations*dt*timedirection)+"\n")
        f.write("iterations\n")
        f.write(str(iterations)+"\n")
        f.write("real runtime (real s)\n")
        f.write(str(time.time()-tstartsecs)+"\n")
        f.write("exit status\n")
        f.write(str(exitstatus(posF)))
        f.close()

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