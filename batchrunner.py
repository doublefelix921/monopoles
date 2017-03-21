from __future__ import print_function, division
import random
import numpy as np
import sys
import time
import trajectory as tj
import usefulfunctions as f
import csv
# import os

# This program runs monopoles.py many times with different
# initial conditions (position, direction of approach, initial 
# kinetic energy). It runs with constant mass, which can be 
# taken in from the command line.

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
settings["StopRadius"] = 20.01 # must start within this radius
settings["StopRadiusSensitivity"] = 0.1  # redo the program with a smaller timestep...
# ...if rho > StopRadius + StopRadiussensitivity.
settings["StopDistance"] = -1 # after traversing an arclength of StopDistance (kpc)...
# ...stop the program even if it hasn't gone past StopRadius. ...
# ...If StopDistance is -1, don't stop based on this variable.
settings["TimeoutSeconds"] = 60
settings["RelativisticVelThreshold"] = 0.15 # stop program if using nonrelativistic regime and  
                                            # the velocity is greater than this * c
settings["PrintStatus"] = True # print info about the program when it starts (later will set to False)

settingsNoPrint = settings # on most runs, don't bother printing the status. This is more commonly used.
settingsNoPrint["PrintStatus"] = False


def convertmany_RectToSpherical(listOfRectangularVecs):
    listOfSphericalVecs = []
    for rectvec in listOfRectangularVecs:
        convertedR = mag(rectvec)
        convertedTheta = np.arccos(rectvec[2]/convertedR)
        convertedPhi = theta_cylindrical(rectvec[0],rectvec[1])
        listOfSphericalVecs += [[convertedR, convertedTheta, convertedPhi]]
    return listOfSphericalVecs

def fibonacci_sphere(samples):
    # this function returns many ~evenly distributed points on a unit sphere.
    # it is taken from stackoverflow: http://stackoverflow.com/a/26127012/3347826
    # I removed the random option in favor of a test that gives a more even distribution.
    # method explained on http://blog.marmakoide.org/?p=1, also saved as pdf in this folder.
    # tests (in the "fibsphere tests" folder) reveal that the points are even to about 8%, 
    # dependent on sample size.

    rnd = 1.
    if samples < 14:
        # tests showed that for samples < 14, this yields a more even distribution.
        rnd = 0.5 * samples
        # rnd = random.random() * samples

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

# RUNNING THE SIMULATIONS
ans = raw_input("Is the data folder empty and ready for batch runs? y/n: ")
if ans!="yes" and ans!="y" and ans!="Yes" and ans!="Y":
    raise Exception("Please make sure the data folder is empty and ready for batch runs.")

# fileList = os.listdir("data")
# if len(fileList)>0:
#     raise Exception("Please make sure the data folder is empty and ready for batch runs.")

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
        except tj.FinalPosError as e:
            print("Must redo run with dt/2 for final position accuracy.")
            print("mag(posF) was %.4gkpc, which exceeds %.4g + %.4g = %.4gkpc." %
                (mag(e.pos), e.settings["StopRadius"], e.settings["StopRadiusSensitivity"],
                    e.settings["StopRadius"] + e.settings["StopRadiusSensitivity"]))
            # print("Starting redo attempt number %d, with a halved timestep.\n" % 
            #     (runcount_finalposerror+1))
            dt = e.dt/2
            runcount_finalposerror += 1
            continue
        except tj.SpeedError as e:
            # redo run with halved timestep
            print("====HAD TO REDO RUN BECAUSE mag(vel) > c ====")
            print("Starting redo attempt number %d with a halved timestep." % 
                (runcount_speederror+1))
            dt = e.dt/2
            runcount_speederror += 1
            continue
    return pos, vel, distanceTraveled, maxVelocity, arcLengthInRho1kpcRegion, iterations

dt = 5.0e10 # timestep, seconds
q_b = 1.066559779e-28 # the dirac magnetic charge of the particle (Ampere*kpc)
m_default = 1e15*1.78266184e-27 # mass (kg, enter in GeV). An estimate from Wick (2002).
try:
    m = float(sys.argv[1])
    print("Mass given as %.4g kg in command line." % m)
except IndexError:
    print("Mass not given in command line. Using default m=%g kg" % m_default)
    m = m_default

#generate positions:
pos0Radius = 20
nPosDirections = 5
points_rectangular = fibonacci_sphere(samples=nPosDirections)
pos0List = [pos0Radius*np.asarray(p) for p in points_rectangular]

#generate kinetic energies:
minKE = (1e11)*1.6827e-49
maxKE = (1e14)*1.6827e-49
nKineticEnergies=2
KE0List = np.linspace(minKE, maxKE, num=nKineticEnergies)

#generate velocities:
nVelDirections = 5

# print status every x runs (starting on 1st)
printInterval = 1

# start with nVelDirections*2 because half point the wrong way
points_rectangular = fibonacci_sphere(samples=nVelDirections*2)
initialConditionList = []
for pos0 in pos0List:
    for KE0 in KE0List:
        vel0DirectionList = [np.asarray(p) for p in points_rectangular if dot(pos0, p) < 0]
        speed = f.speed_from_KE(KE0, m)
        for vel0Dir in vel0DirectionList:
            initialConditionList += [[pos0, speed*vel0Dir, dt, m, q_b, settingsNoPrint]]

outputData = []
lenICList = len(initialConditionList)
for n in range(lenICList):
    c = initialConditionList[n]
    if n%printInterval == 0:
        print("\nBeginning run %d of %d..." % (n+1, lenICList), end="")
        c[-1] = settings
    outputData += [run_trajectory_wrapper(c[0], c[1], c[2], c[3], c[4], c[5])]
    if n%printInterval ==0: 
        print("done.")

# write to a CSV, yo.
timestamp = time.strftime("%Y.%m.%d.%H.%M.%S")
csvrows = []
csvrows += [["SETTINGS:"]]
csvrows += [["SimulateBool", "DiskPlotBool", "TrajectoryBool", "VecPlotBool", "XYCrossBool",
    "XZCrossBool", "BFieldTests", "TimeDirection", "Relativistic", "HaloFieldBool", "XFieldBool", 
    "DiskFieldBool", "StopRadius", "StopRadiusSensitivity", "StopDistance", "TimeoutSeconds", 
    "RelativisticVelThreshold", "PrintStatus"]]
settingsRowIndex = csvrows.index(["SETTINGS:"])
csvrows += [[str(settings[setting]) for setting in csvrows[settingsRowIndex+1]]]
csvrows += [[]]

# NOW ADD INITIAL AND FINAL CONDITIONS AND RUN NUMBER FOR EACH RUN
csvrows += [["RUN DATA:"]]
csvrows += [["run #", "mass", "time (Y.m.d.H.M.S)", "dt", "q_b", 
    "pos0[x]", "pos0[y]", "pos0[z]", "vel0[x]", "vel0[y]", "vel0[z]",
    "posF[x]", "posF[y]", "posF[z]", "velF[x]", "velF[y]", "velF[z]", "distanceTraveled", 
    "maxVelocity", "arcLengthInRho1kpcRegion", "iterations"] ]
print("Formatting data to write to a .csv... ",end="")
for n in range(len(initialConditionList)):
    IC = initialConditionList[n]
    FC = outputData[n]
    csvrows += [[n+1, m, timestamp, dt, q_b, IC[0][0], IC[0][1], IC[0][2], IC[1][0], IC[1][1], IC[1][2],
    FC[0][0], FC[0][1], FC[0][2], FC[1][0], FC[1][1], FC[1][2], FC[2], FC[3], FC[4], FC[5] ] ]
print("done.")

outputfilename = "data/batchrun" + timestamp + ".csv"
print("Writing to file %s (labeled by time of first run)... " % outputfilename, end="")
with open(outputfilename, 'wb') as file:
    mywriter = csv.writer(file)
    for line in csvrows:
        mywriter.writerow(line)
    file.close()
print("done.\n")

# initialConditionList has pos0, vel0, dt, m, qb, settings
# outputdata has pos, vel, distanceTraveled, maxVelocity, arcLengthInRho1kpcRegion, iterations

# DurationMins = (time.time() - startTime)/60
# print("\nDone with %d of %d runs." % (numRuns, numRuns))
# print("Took %.4g minutes, or %.5gs per run." % (DurationMins, DurationMins*60/numRuns))

