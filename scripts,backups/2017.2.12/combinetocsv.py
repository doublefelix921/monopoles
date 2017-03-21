from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass
try: range = xrange
except NameError: pass

import os
import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import datetime
import csv

# what this program does:
# 1. combine data files of completed runs in the data/ folder to a single csv
# 2. check the min/max of final rho values
# 3. add failed runs to the csv; it will also try to add an error code in 
#    the 1st row (index 0) and will be printed after this program completes.
#    Error codes:
#       "POS0ERROR" - the initial position was too far from the center and 
#       it triggered the exit feature

def txtToNpArray(text):
    text = text.replace("[","")
    text = text.replace("]","")
    # text = text.replace(" ","")
    text = text.split()
    text2=[]
    for n in text:
        # print("converting to float: %s"%n)
        text2+=[float(n)]
    return np.asarray(text2)
 
def mag(V):
    return np.sqrt(V[0]**2+V[1]**2+V[2]**2)

# EXTRACTING THE DATA FROM FILES
filelist = os.listdir("data")
answer = raw_input("Is the data folder ready for ALL files to be combined " \
    "into a single csv? \nFiles will not be deleted. \ny/n: ")
if answer == "no" or answer == "n" or answer == "N":
    print("exiting.")
    sys.exit()

finalpoints = [] #contains final E, theta, phi
finalrhos = []

unfinishedfiles = []
failedruns = []

outputdata = [] #contains everything for all runs

for n in range(len(filelist)):
    print("reading data/%s..." % filelist[n])
    f = open("data/"+filelist[n], "r+")
    lines = f.readlines()
    lines = [x.replace("\n","") for x in lines]
    row = []

    # COL DESCRIPTIONS
    # this is also the order in which variables are written in a 
    # singlerun-file but they do not have an error code except in the case
    # of a POS0ERROR, in which case that code is at the end of the 
    # (incomplete) file.
    colheaders = ["ERROR CODE", "start time (YYYYMMDD-hhmmssms)", "v0 (kpc/s)", 
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

    #INITIAL CONDITIONS
    # str, str, nparray, float, float, float, float, float, float, float, 
    # float, float
    row += [lines[2]]   #start time (str)
    row += [lines[4]]   # v0
    row += [lines[6]]   
    row += [lines[8]]   
    row += [lines[10]]  # E0
    row += [lines[12]] 
    row += [lines[14]] 
    row += [lines[16]] 
    row += [lines[18]] 
    row += [lines[20]] 
    row += [lines[22]] # distance to halt (kpc)

    if lines[23][:5]=="ERROR":
        if lines[24] == "POS0ERROR":
            row[0] = "POS0ERROR"
            failedruns+=[row]
            continue

    # FINAL CONDITIONS
    # str, nparray, float, float, nparray, float, float, float, float, 
    # float, float, float, float, float, float, int, float, nparray, bool
    try:
        row += [lines[25]]  # finish time
        row += [lines[27]]
        row += [lines[29]]  # mag(vel)
        row += [lines[31]]
        row += [lines[33]]  # pos
        finalrhos += [mag(txtToNpArray(lines[33]))]
        row += [lines[35]]
        row += [lines[37]]
    except IndexError:
        print("Skipping file %s as it never finished." % filelist[n])
        unfinishedfiles += [filelist[n]]
        continue

    
    row += [lines[39]]
    row += [lines[41]]
    row += [lines[43]]
    row += [lines[45]]
    row += [lines[47]]
    row += [lines[49]]
    row += [lines[51]]
    row += [lines[53]]
    row += [lines[55]]
    row += [lines[57]]
    row += [lines[59]]  #final acceleration
    row += [lines[63]]  #exit status
    row = ["SUCCESS"] + row
    print("adding this to outputdata:",row)
    outputdata+=[row]
print("currently, outputdata is", outputdata)
outputdata = [colheaders] + outputdata
print("done.")
print("First run in the batch started at time %s (YYYYMMDD-hhmmssms)." % 
    outputdata[1][1])

outputfilename = "data/batchrun" + outputdata[1][1] + ".csv"
print("Writing to file %s (labeled by time of first run)... " % outputfilename, end="")

with open(outputfilename, 'wb') as csvfile:
    mywriter = csv.writer(csvfile)
    for line in outputdata:
        mywriter.writerow(line)
    csvfile.close()
print("done.\n")

print("===STATISTICS FOR ALL %d RUNS===" % (len(outputdata)-1))
print("Final 'rho' values ranged from %.2g to %.2g" % 
    (min(finalrhos),max(finalrhos)))

print("%d files only had start data." % len(unfinishedfiles), end="")
if len(unfinishedfiles) > 0:
    print(" Their filenames are:")
else: print()
for n in range(len(unfinishedfiles)):
    print("    "+unfinishedfiles[n])

errors = [[n,outputdata[n][0]] for n in range(len(outputdata)) if 
    (outputdata[n][0]!="SUCCESS" and outputdata[n][0]!="ERROR CODE")]


if len(errors)>0:
    print("%d runs had errors. They are:" % len(errors))
    print("    Run # | outputdata[n] | error code")
    for n in range(len(errors)):
        print("    %-8d%-16d%s" % (errors[n][0]+1, errors[n][0], errors[n][1]))
else:
    print("0 runs had errors.")
    
# PLOTTING THE FINAL POSITIONS
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# xs = [p[0] for p in finalpoints]
# ys = [p[1] for p in finalpoints]
# zs = [p[2] for p in finalpoints]
# ax.scatter(xs, ys, zs)

# ax.set_xlabel('X')   #label each positional axis
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')

# plt.show()


# f.close()