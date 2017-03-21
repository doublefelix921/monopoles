from __future__ import print_function, division
import math, random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import os
import time

# This program is used to run monopoles.py many times with different
# initial conditions (position, direction of approach, initial 
# kinetic energy). monopoles.py takes the following cmdline args:
    # E0    = sys.argv[1]  #E0 is initial kinetic energy
    # vtheta0  = sys.argv[2] #vtheta0 and vphi0 give the monopole's initial 
    # vphi0    = sys.argv[3] #direction of approach
    # x0      = sys.argv[4]
    # y0      = sys.argv[5]
    # z0      = sys.argv[6]

# this program should take in mass from the command line.
# for starters, let's use very un-fine changes in variables.

combineToCSV = True

def mag(V):
    return np.sqrt(V[0]**2+V[1]**2+V[2]**2)

def dot(u, v):
    #u and v must be in cartesian coordinates
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

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

def convertmany_RectToSpherical(listOfRectangularVecs):
    listOfSphericalVecs = []
    for rectvec in listOfRectangularVecs:
        convertedR = mag(rectvec)
        convertedTheta = np.arccos(rectvec[2]/convertedR)
        convertedPhi = Theta(rectvec[0],rectvec[1])
        listOfSphericalVecs += [[convertedR, convertedTheta, convertedPhi]]
    return listOfSphericalVecs

def fibonacci_sphere(samples=1,randomize=True):
    # this function is taken from stack exchange.
    # credit should be given to them for it, not me.
    # their content, I believe, goes under the MIT license.
    # returns many ~evenly (randomly?) distributed points on a unit sphere
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points_rectangular = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points_rectangular.append([x,y,z])

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
    print("exiting. Please ready the data folder.")
    sys.exit()
filelist = os.listdir("data")
if len(filelist)>0:
    print("Data folder is not empty. " \
        "Please ready the data folder and try again.")
    sys.exit(1)

# convert to spherical coordinates
# positions_spherical = []
# for p in points_rectangular:
#     theta = np.arccos(p[2]/mag(p))
#     phi = Theta(p[0], p[1])
#     positions_spherical += [[1, theta, phi]]
# print("created %d points." % len(positions_spherical))


# CALL A RUN FROM A VARIETY OF DIRECTIONS, STARTING AT A SPECIFIED PT:
def runConstantPos0_E0_andVariableDirection(pt0, E0, nDirections):
    #possible directions are only those which would enter the sphere.
    #thus the dot product with the initial position pt0 should be negative.
    points_rectangular = fibonacci_sphere(samples=nDirections*2)

    dirs = [p for p in points_rectangular if dot(pt0,p) < 0]

    # PLOTTING THE DIRECTIONS
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # xs = [p[0] for p in dirs]
    # ys = [p[1] for p in dirs]
    # zs = [p[2] for p in dirs]
    # ax.scatter(xs, ys, zs)
    # ax.set_xlabel('X (kpc)')
    # ax.set_ylabel('Y (kpc)')
    # ax.set_zlabel('Z (kpc)')
    # plt.show()

    dirs = convertmany_RectToSpherical(dirs)
    numruns = len(dirs)
    starttime = time.time()
    for n in range(numruns):
        print("\nStarting run %d of %d..." % (n+1, numruns))
        runstring = "python monopoles.py " + str(E0)+" " + str(dirs[n][1]) + " " + str(dirs[n][2]) + " "
        runstring += str(pt0[0]) + " " + str(pt0[1]) + " " + str(pt0[2]) 
        # print("\nrunning with command:\n'%s'\n" % runstring)
        os.system(runstring)
    print("\nDone with %d of %d runs." % (numruns, numruns))
    durationmins = (time.time() - starttime)/60
    print("Took %g minutes." % durationmins)

runConstantPos0_E0_andVariableDirection([20, 0, 0], (1e13)*1.6827e-49,10)

# def runConstantE0_andVariablePos0_Direction(E0, nPositions, nDirections):
    # runs from all position + direction combinations.
    # space of directions is half as large as that of positions, so I think it's
    # advisable to let nDirections ~= (1/2)*nPositions.

if combineToCSV:
    with open("combinetocsv.py") as f:
        code = compile(f.read(), "combinetocsv.py", "exec")
        exec(code)

# runs starting from earth, varying all components of v0.
# (mapping flux out to flux on sphere)
# for ps in positions_spherical:
#     E0 = (1e13)*1.6827e-49 #enter in GeV
#     # runstring = "python monopoles.py %f %f %f " % (E0, p[1], p[2])
#     runstring = "python monopoles.py "+str(E0)+" "+str(ps[1])+" "+str(ps[2])
#     print("running with command:\n'%s'" % runstring)
#     os.system(runstring)

