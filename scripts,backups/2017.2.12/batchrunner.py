from __future__ import print_function, division
import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import os
import time
import trajectory as tj
import usefulfunctions as f

# This program is used to run monopoles.py many times with different
# initial conditions (position, direction of approach, initial 
# kinetic energy). 

# this program should take in mass from the command line.
# for starters, let's use very un-fine changes in variables.

combineToCSV = True

mag = f.mag
dot = f.dot
theta_cylindrical = f.theta_cylindrical

settings = {}
settings["SimulateBool"] = True    # determines if trajectory simulation should occur.
settings["DiskPlotBool"] = False  # determines if the various plots are made.
settings["TrajectoryBool"] = False
settings["VecPlotBool"] = False
settings["XYCrossBool"] = False
settings["XZCrossBool"] = False
settings["BFieldTests"] = False # if true, find avg & min/max of bfield.
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




def convertmany_RectToSpherical(listOfRectangularVecs):
    listOfSphericalVecs = []
    for rectvec in listOfRectangularVecs:
        convertedR = mag(rectvec)
        convertedTheta = np.arccos(rectvec[2]/convertedR)
        convertedPhi = theta_cylindrical(rectvec[0],rectvec[1])
        listOfSphericalVecs += [[convertedR, convertedTheta, convertedPhi]]
    return listOfSphericalVecs

def fibonacci_sphere(samples=1, randomize=True):
    # this function is taken from stack exchange.
    # credit should be given to them for it, not me.
    # their content, I believe, goes under the MIT license.
    # returns many ~evenly (randomly?) distributed points on a unit sphere
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points_rectangular = []
    offset = 2./samples
    increment = np.pi * (3. - np.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = np.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = np.cos(phi) * r
        z = np.sin(phi) * r

        points_rectangular.append([x, y, z])

    return points_rectangular

# PLOTTING THE DIRECTIONS
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# xs = [p[0] for p in points_rectangular]
# ys = [p[1] for p in points_rectangular]
# zs = [p[2] for p in points_rectangular]
# ax.scatter(xs, ys, zs)
# ax.set_xlabel('X (kpc)')
# ax.set_ylabel('Y (kpc)')
# ax.set_zlabel('Z (kpc)')
# plt.show()

# RUNNING THE SIMULATIONS
ans = raw_input("Is the data folder empty and ready for batch runs? y/n: ")
if ans!="yes" and ans!="y" and ans!="Yes" and ans!="Y":
    raise Exception("Please make sure the data folder is empty and ready for batch runs.")

fileList = os.listdir("data")
if len(fileList)>0:
    raise Exception("Please make sure the data folder is empty and ready for batch runs.")
# convert to spherical coordinates
# positions_spherical = []
# for p in points_rectangular:
#     theta = np.arccos(p[2]/mag(p))
#     phi = theta_cylindrical(p[0], p[1])
#     positions_spherical += [[1, theta, phi]]
# print("created %d points." % len(positions_spherical))

def run_trajectory_wrapper(pos0, vel0, dt, m, q_b, settings):
    runcount_finalposerror = 0
    runcount_speederror = 0
    successfulRun = False
    while successfulRun == False:
        try:
            pos, vel, distanceTraveled, maxVelocity, \
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
    return pos, vel, distanceTraveled, maxVelocity, arcLengthInRho1kpcRegion, iterations


# CALL A RUN FROM A VARIETY OF DIRECTIONS, STARTING AT A SPECIFIED PT:
def run_variable_velhat(pos0, KE0, nVelDirections, m, dt, q_b, settings):
    #possible directions are only those which would enter the sphere.
    #thus the dot product with the initial position pt0 should be negative.
    points_rectangular = fibonacci_sphere(samples=nVelDirections*2)
    velocityDirections = [np.asarray(p) for p in points_rectangular if dot(pos0, p) < 0]

    speed = f.speed_from_KE(KE0, m)
    outcomeData_variable_velhat = []
    numRuns = len(velocityDirections)
    startTime = time.time()
    for n in range(numRuns):
        print("Run %d of %d: " % (n+1, numRuns), end="")
        outcomeData_variable_velhat += [[
        run_trajectory_wrapper(pos0, speed*velocityDirections[n], dt, m, q_b, settings)]]

    DurationMins = (time.time() - startTime)/60
    print("\nDone with %d of %d runs." % (numRuns, numRuns))
    print("Took %.4g minutes, or %.5gs per run." % (DurationMins, DurationMins*60/numRuns))
    return outcomeData_variable_velhat

def run_variable_poshat_velhat(nPosDirections, KE0, nVelDirections, m, dt, q_b, settings):
    # generate the pos directions:
    posDirections = fibonacci_sphere(samples=nPosDirections)
    outcomeData_variable_poshat_velhat = [
        run_variable_velhat(pos, KE0, nVelDirections, m, dt, q_b, settings) for pos in posDirections]
    return outcomeData_variable_poshat_velhat

# outcomeData_variable_velhat = run_variable_velhat(pos0=[20., 0., 0.], 
#     KE0=(1e13)*1.6827e-49, nVelDirections=10, 
#     m=1e15*1.78266184e-27, dt=5e10, q_b=1.06655979e-28, settings=settings)

outcomeData_variable_poshat_velhat = run_variable_poshat_velhat(nPosDirections=10, 
    KE0=(1e13)*1.6827e-49, nVelDirections=10, 
    m=1e15*1.78266184e-27, dt=5e10, q_b=1.06655979e-28, settings=settings)

if 

# if combineToCSV:
#     with open("combinetocsv.py") as f:
#         code = compile(f.read(), "combinetocsv.py", "exec")
#         exec(code)

# runs starting from earth, varying all components of v0.
# (mapping flux out to flux on sphere)
# for ps in positions_spherical:
#     E0 = (1e13)*1.6827e-49 #enter in GeV
#     # runstring = "python monopoles.py %f %f %f " % (E0, p[1], p[2])
#     runstring = "python monopoles.py "+str(E0)+" "+str(ps[1])+" "+str(ps[2])
#     print("running with command:\n'%s'" % runstring)
#     os.system(runstring)

