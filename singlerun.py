from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass
try: range = xrange
except NameError: pass

import time
import numpy as np
import trajectory as tj
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


settings = {}
settings["SimulateBool"] = True    # determines if trajectory simulation should occur.
settings["DiskPlotBool"] = True  # determines if the various plots are made.
settings["TrajectoryBool"] = True
settings["VecPlotBool"] = True
settings["XYCrossBool"] = False
settings["XZCrossBool"] = False
settings["BFieldTests"] = True # if true, find avg & min/max of bfield.
settings["TimeDirection"] = 1   # track forwards in time if 1, backwards if -1.
settings["Relativistic"] = True # decides whether to use Relativistic or nonRelativistic...
# ...expression for acceleration. nonRelativistic computes much more quickly.
settings["HaloFieldBool"] = True     # turn on or off various components of the field.
settings["XFieldBool"] = True        # the random field is not used as it averages out
settings["DiskFieldBool"] = True          # anyway for this application.
settings["StopRadius"] = 20.01 # stop the program if rho > StopRadius.
settings["StopRadiusSensitivity"] = 0.1  # redo the program with a smaller timestep...
# ...if rho > StopRadius + StopRadiussensitivity.
settings["StopDistance"] = -1 # after traversing an arclength of StopDistance (kpc)...
# ...stop the program even if it hasn't gone past StopRadius. ...
# ...If StopDistance is -1, don't stop based on this variable.
settings["TimeoutSeconds"] = 60
settings["RelativisticVelThreshold"] = 0.15 # stop program if using nonrelativistic regime and  
                                            # the velocity is greater than this * c
settings["PrintStatus"] = True

# =======CONSTANTS=======
c = tj.c # speed of light (seconds)
kgkpc2s2_To_GeV = tj.kgkpc2s2_To_GeV

# =======DEFAULT STARTING CONDITIONS=======
q_b = 1.066559779e-28 # the dirac magnetic charge of the particle (Ampere*kpc)
m = 1e15*1.78266184e-27 # mass (kg, enter in GeV). An estimate from Wick (2002).
pos_earth = np.array([-8.5, 0, 0]) # position of earth is x=-8.33937979kpc but we must...
# ...use -8.5 because the model is based on it.
# pos0 = np.array([-19.8, 1.5, 1])
pos0 = pos_earth # starting position (kpc), 
# v_toEarth = np.array([-8.5,0,0]) - pos0       # point velocity to the earth
# vel0 = c*(.4e-3)*(v_toEarth)*(1.0/mag(v_toEarth))
E0 = (10**11)*1.6827e-49 # E0, enter in GeV in the parenthesis
vel0mag = E0*c/np.sqrt(c**4*m**2 + E0**2)
vtheta0 = 1.7
vphi0 = 3.14
vel0 = vel0mag*f.sph_to_rect(1, vtheta0, vphi0)
dt_default = 5.0e10
dt = dt_default

print("Note: Run me with python -i singlerun.py to keep terminal " \
    "open afterwards.")

if not settings["HaloFieldBool"]:
    print("WARNING: The halo field component is disabled (halobool==False)")
if not settings["XFieldBool"]:
    print("WARNING: The x field component is disabled (xfieldbool==False)")
if not settings["DiskFieldBool"]:
    print("WARNING: The disk field component is disabled (diskbool==False)")
if not settings["SimulateBool"]:
    if settings["TrajectoryBool"]:
        print('ERROR: settings["TrajectoryBool"]==True but settings["SimulateBool"]==False ! Halting.\n')
        sys.exit(0)
    print('WARNING: You have chosen not to simulate the trajectory. ' +
        '(settings["SimulateBool"]==False)')

if settings["SimulateBool"]:
    runcount_finalposerror = 0
    runcount_speederror = 0
    SuccessfulRun = False
    while SuccessfulRun == False:
        try:
            posF, velF, distanceTraveled, maxVelocity, \
                arcLengthInRho1kpcRegion, iterations = \
                tj.run_trajectory(pos0, vel0, dt, m, q_b, settings)
            successfulRun = True
        except tj.FinalPosError:
            print("ERROR: Must redo run for final position accuracy.")
            print("mag(posF) was %.4gkpc, which is greater than the " \
                "threshhold of %.4g + %.4g = %.4gkpc." %
                (tj.mag(posF), settings["StopRadius"], settings["StopRadiusSensitivity"],
                    settings["StopRadius"]+settings["StopRadiusSensitivity"]))
            print("Starting redo attempt number %d, with a halved timestep.\n" % 
                (runcount_finalposerror+1))
            dt = dt/2
            runcount_finalposerror += 1
            continue
        except tj.SpeedError:
            # redo run with halved timestep
            print("====HAD TO REDO RUN BECAUSE mag(vel) > c ====")
            print("Starting redo attempt number %d with a halved timestep." % 
                (runcount_speederror+1))
            dt = dt/2
            runcount_speederror += 1
            continue

    # if you are up to here, the run has succeeded.
    reruns = runcount_finalposerror + runcount_speederror

    # print some stuff after the run?
    print("Run complete with no errors after %d reruns." % reruns)
    print("iterations:", iterations)
    print("posF:   ", posF)
    print("velF/c: ", velF/c)
    print("Initial KE(GeV): %.4g" % (tj.KE(m, vel0) * kgkpc2s2_To_GeV))
    print("Final KE(GeV): %.4g" % (tj.KE(m, velF) * kgkpc2s2_To_GeV))

if settings["BFieldTests"]==True:
    with open("bfieldtests.py") as f:
        code = compile(f.read(), "bfieldtests.py", "exec")
        print()
        exec(code)

if (settings["DiskPlotBool"] or \
    settings["TrajectoryBool"] or \
    settings["XYCrossBool"] or \
    settings["XZCrossBool"]) \
    and not settings["TrajectoryBool"]:
    # trajectory was already done in the run_trajectory function,
    # as it needed local variables. All other plots were also done at the time.
    # Only do this if they haven't been done yet.
    with open("grapher.py") as f:
        code = compile(f.read(), "grapher.py", 'exec')
        print()
        exec(code)
print()