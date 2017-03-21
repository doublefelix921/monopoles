from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass
try: range = xrange
except NameError: pass

import time
import numpy as np
import sys
import random
import getB as gb
import usefulfunctions as f

# Tracking Magnetic Monopoles Through the Galaxy:
# input starting conditions into the program and it will track the magnetic
# monopole's path through the Bfield given by "A New Model of the Galatic 
# Magnetic Field" Jansson, Farrar (2012).

# units used are all SI units except: distance is in kpc, and angles are in
# degrees, magnetic field strength inputs are in microGauss but calculations
# are done in Tesla so the gb.get_B() function does conversions before output.

# variable names generally match those in the Bfield paper.

# physics variables are in their traditional case (usually lower).
# ^ this has prescedence over other naming conventions.
# settings variables are in CamelCaseLikeThis
# functions should be lower_case_with_underscore_separations

# =======CONSTANT PARAMETERS=======
c = f.c            # light speed in kpc/s

# a dict with the run's settings. Initialize for no reason now.
settingsNames = ["SimulateBool", "DiskPlotBool", 
"TrajectoryBool", "VecPlotBool", "XYCrossBool", "XZCrossBool", "BFieldTests",
"TimeDirection", "Relativistic", "HaloFieldBool", "XFieldBool", "DiskFieldBool",
"StopRadius", "StopRadiusSensitivity", "StopDistance", "TimeoutSeconds",
"RelativisticVelThreshold", "PrintStatus"]

# magnitude of a vector
mag = f.mag

# rectangular coordinates (x, y) to cylindrical (azimuthal) angle
theta_cylindrical = f.theta_cylindrical

def exit_status(pos, stopRadius):
    if mag(pos) >= stopRadius:
        return True

# The lorentz factor
def gamma(vel):
    magvel = mag(vel)
    return c/(np.sqrt(c+magvel)*np.sqrt(c-magvel))

# Kinetic Energy
def KE(m, vel):
    #return in kg*kpc^2/s^2
    if gamma(vel) > 1.000000001:
        return m*c**2*(gamma(vel)-1)
    else:
        return m*mag(vel)**2/2

#returns Relativistic acceleration. only necessary for large velocity.
def acc_relativistic(pos, vel, q_b, m):   
    bfield = gb.get_B(pos)
    force = q_b*bfield    # In kg*kpc/s^2
    acc = (1.0/(gamma(vel)*m*c**2)) * (c**2*force - np.dot(force, vel)*vel)
    return acc

def vel(pos, vel):
    return vel

#this is really only a function of position.
def acc_classical(pos, vel, q_b, m):
    bfield = gb.get_B(pos)
    acc  = bfield*(q_b/m)
    return acc

def euler_step(pos, vel, drdt, dvdt, dt):
    dr = drdt(pos,vel)*dt
    dv = dvdt(pos,vel)*dt
    return dr, dv
    
def rk4_step(pos, vel, drdt, dvdt, dt, acc_function, q_b, m):
    #return the step to be added, NOT the updated vector, so that dr can be 
    #added to distance_traveled.
    k1 = np.array([vel*dt,                    acc_function(pos,           vel, q_b, m)*dt           ])
    k2 = np.array([(vel+k1[1]/2.0)*dt,        acc_function(pos+k1[0]/2.0, vel+k1[1]/2.0, q_b, m)*dt ])
    k3 = np.array([(vel+k2[1]/2.0)*dt,        acc_function(pos+k2[0]/2.0, vel+k2[1]/2.0, q_b, m)*dt ])
    k4 = np.array([(vel+k3[1])*dt,            acc_function(pos+k3[0],     vel+k3[1], q_b, m)*dt     ])
    
    pos_step = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
    # tested: mutability of np array does not mess this up (changing vel, pos does not chg k1)
    vel_step = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0
    
    return pos_step, vel_step

class Pos0Error(Exception):
    def __init__(self, rho, stopRadius, stopRadiusSensitivity, iterations):
        self.message = "Pos0Error on iteration %d: Initial position was outside the sim volume!\n" + \
            "Rho was %.4g and sim radius is" % (iterations, rho, stopRadius + stopRadiusSensitivity)
    def __str__(self):
        return repr(self.message)

class SpeedError(Exception):
    def __init__(self, vel, iterations, dt):
        self.vel = vel
        self.iterations = iterations
        self.dt = dt
        self.message = "SpeedError on iteration %d: mag(vel) >= c !  " + \
            "speed/c = %.4g,   gamma = %.4g" % (iterations, mag(vel)/c, gamma(vel))
    def __str__(self):
        return repr(self.message)

class TimeoutError(Exception):
    def __init__(self, TimeElapsedSeconds, TimeoutSeconds, iterations):
        self.message = "TimeoutError on iteration %d: " + \
        "One run lasted %ds, which is more than the allowed time of " + \
        "settings['TimeoutSeconds']= %ds \n" % (iterations, TimeElapsedSeconds, TimeoutSeconds)
    def __str__(self):
        return repr(self.message)

class FinalPosError(Exception):
    def __init__(self, pos, stopRadius, stopRadiusSensitivity, iterations, settings, dt):
        self.pos = pos
        self.stopRadius = stopRadius
        self.stopRadiusSensitivity = stopRadiusSensitivity
        self.iterations = iterations
        self.settings = settings
        self.dt = dt
        stringg = ("FinalPosError on iteration %d: Must redo run for final position accuracy.\n" + \
            "mag(posF) was %.4gkpc, which is greater than the " + \
            "threshold of %.4g + %.4g = %.4gkpc.") % \
            (iterations, mag(pos), stopRadius, stopRadiusSensitivity, (stopRadius + stopRadiusSensitivity))
        self.message = stringg
    def __str__(self):
        return repr(self.message)

class RelativisticRegimeError(Exception):
    def __init__(self, vel, Threshold, iterations):
        self.message = "RelativisticRegimeError on iteration %d: Was using acc_classical, " + \
        "but mag(vel)/c is %.4g, and gamma is %.4g! " + \
        "This is beyond the threshold of %.4g c." %  (iterations, mag(vel)/c, gamma(vel), threshold)
    def __str__(self):
        return repr(self.message)

# idea: make a Relativistic and nonrel paradigm. when the speed is nonrel (difference
# would be less than machine precision) Relativistic=False. Every X iterations, check if 
# the paradigm should be changed, where X is the least number of iterations it would take
# to change the velocity by, say, 1% of the speed of light.

# =======PROGRAM BEGINS=======
def run_trajectory(pos0, vel0, dt, m, q_b, settings):
    # ======== PART 1: PREPARATION ========    
    #==optional parameters & other stuff==
    tStartSecs = time.time()
    maxVelocity = mag(vel0)
    iterations = 0
    arcLengthInRho1kpcRegion = 0.0
    
    if exit_status(pos0, settings["StopRadius"]) == True:
        raise Pos0Error(mag(pos0), settings["StopRadius"], 
            settings["StopRadiusSensitivity"], iterations)

    #for graphing
    xList, yList, zList = [np.array([]) for _ in range(3)]

    if settings["Relativistic"]: acc_function = acc_relativistic
    else:            acc_function = acc_classical

    if settings["HaloFieldBool"]: gb.HaloFieldBool=True
    if settings["XFieldBool"]:    gb.XFieldBool   =True     # the random field is not used as it averages out
    if settings["DiskFieldBool"]: gb.DiskFieldBool=True     # anyway for this application.

    # ======== PART 2: LOOP ========
    pos = np.copy(pos0)
    vel = np.copy(vel0)
    acc = acc_function(pos, vel, q_b, m) # python complains unless acc is defined 
                                # outside the loop.
    distanceTraveled = 0.0   #halt the loop if the monopole travels too far.

    if settings["PrintStatus"]  == True:
        print("Starting traversal. Time is moving ",end="")
        if settings["TimeDirection"]==1:  print("FORWARDS.")
        if settings["TimeDirection"]==-1: print("BACKWARDS.")
        if settings["Relativistic"] == True:  print("Using Relativistic acceleration.")
        if settings["Relativistic"] == False: print("Using CLASSICAL acceleration.")

    while mag(pos) < settings["StopRadius"]:
        iterations += 1

        # update trajectory. Use settings["TimeDirection"] =-1 to track backwards.
        acc = acc_function(pos, vel, q_b, m)
        pos_step, vel_step = rk4_step(pos, vel, vel, acc, dt, acc_function, q_b, m)
        a = settings["TimeDirection"]*pos_step
        try: 
            pos += settings["TimeDirection"]*a
        except:
            print("settings timedirection:", settings["TimeDirection"], type(settings["TimeDirection"]))
            print("pos_step:", pos_step, type(pos_step))
            print("a:",a,type(a))
            print("pos", pos, type(pos))
            sys.exit(1)
        vel += settings["TimeDirection"]*vel_step
        distanceTraveled += mag(pos_step)
        
        if settings["TrajectoryBool"]:
            # FOR GRAPHING:
            xList = np.concatenate((xList, np.array([pos[0]])))
            yList = np.concatenate((yList, np.array([pos[1]])))
            zList = np.concatenate((zList, np.array([pos[2]])))

        if mag(vel) >= c: raise SpeedError(vel, iterations, dt)
        if iterations % 10000 == 0:
            if time.time()-tStartSecs > settings["TimeoutSeconds"]:
                raise TimeoutError(time.time() - tStartSecs, settings["TimeoutSeconds"], iterations)
            bf = gb.get_B(pos)
            disp = pos - pos0
            print("============done with iteration no.", iterations, "============")
            print("pos-pos0 = [%.4g, %.4g, %.4g]. |pos-pos0| = %.4gkpc" % (disp[0], disp[1], disp[2],mag(disp)))
            print("vel/c = [%.4g, %.4g, %.4g]. |vel|/c = %.4g" % (vel[0]/c, vel[1]/c, vel[2]/c,mag(vel)/c))            
            print("acc/c = [%.4g, %.4g, %.4g]. |acc|/c = %.4g/s" % (acc[0]/c, acc[1]/c, acc[2]/c,mag(acc)/c))            
            print("bfield = [%.4g, %.4g, %.4g]. |bfield| = %.4gT" % (bf[0], bf[1], bf[2],mag(bf)))
            print("arc length = %.4gkpc" % distanceTraveled)
            print("simulated time = %.4gs" % (iterations*settings["TimeDirection"]*dt))
            print("real runtime = %.4g real seconds" % (time.time()-tStartSecs))
            print()
            if settings["Relativistic"] == False and mag(vel) > 0.2*c:
                print("\n=====================================================")
                print("====================ERROR===========================")
                print("=== mag(velocity) > %.2fc but you are not using    ===" % \
                    settings["RelativisticVelThreshold"])
                print("===      the Relativistic acceleration!!!         ===")
                raise RelativisticRegimeError(vel, settings["RelativisticVelThreshold"], iterations)
            
        if mag(pos) < 1:
            arcLengthInRho1kpcRegion += mag(pos_step)
    if settings["PrintStatus"] == True:
        if exit_status(pos, settings["StopRadius"]) == True:
            print("Stopped because mag(pos)=%.4g >= %.4g.\n" % 
                (mag(pos), settings["StopRadius"]))
        if (distanceTraveled > settings["StopDistance"]) and settings["StopDistance"] > 0:
            print("Stopped because distance traveled > %.4g. \n" % settings["StopDistance"])

    if mag(pos) > settings["StopRadius"] + settings["StopRadiusSensitivity"]:
        raise FinalPosError(pos, settings["StopRadius"], 
            settings["StopRadiusSensitivity"], iterations, settings, dt)
    distance_from_start = mag(pos-pos0)
    thetaF = np.arccos(pos[2]/mag(pos))
    phiF = theta_cylindrical(pos[0], pos[1])

    if settings["TrajectoryBool"]:
        with open("grapher.py") as f:
            code = compile(f.read(), "grapher.py", 'exec')
            print()
            exec(code)

    return pos,vel, distanceTraveled, maxVelocity, \
        arcLengthInRho1kpcRegion, iterations

print()