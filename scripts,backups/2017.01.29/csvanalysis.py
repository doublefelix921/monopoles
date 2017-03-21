from __future__ import print_function, division
try: input = raw_input
except NameError: pass

import csv
import sys
import matplotlib.pyplot as plt
import numpy as np

c = 9.715611713408621e-12

# RUNDATA COLUMN HEADERS
# [0]  ERROR CODE
# [1]  start time (YYYYMMDD-hhmmssms)
# [2]  v0 (kpc/s)
# [3]  mag(v0) (kpc/s)
# [4]  mag(v0/c)
# [5]  E0 (kg*kpc^2/s^2)
# [6]  theta0 (rads)
# [7]  phi0 (rads)
# [8]  q_b (A*kpc)
# [9]  mass (kg)
# [10] dt (s)
# [11] distance to halt (kpc)
# [12] end time (YYYYMMDD-hhmmssms)
# [13] vel_final (kpc/s)
# [14] mag(vel_final) (kpc/s)
# [15] mag(vel_final/c)
# [16] pos_final(kpc)
# [17] mag(pos_final) (kpc)
# [18] theta_final (rads)
# [19] phi_final (rads)
# [20] distance tracked (kpc)
# [21] distance from start (kpc)
# [22] KE_final (kg*kpc^2/s^2)
# [23] max velocity (kpc/s)
# [24] max velocity/c
# [25] arclength in rho1kpc region (kpc)
# [26] time after t=0 start (s)
# [27] iterations
# [28] real runtime (real s)
# [29] acc_final (kpc/s/s)
# [30] exit status (True/False)


#CSV importer.
#does not handle commas or newlines within entries.
#scroll to the bottom to enter CSV files to import & clean

def transpose(data):
    trans = [list(x) for x in zip(*data)]
    return trans

def convertToNestedList(data):
    print("==> Converting data to nested list... ", end="")
    newdata = []
    for n in range(len(data)):
        newdata+=[data[n].split(",")]
    print("done.")
    print("Data has %d rows. Its first, second, and last lines " \
        "have %d, %d, and %d columns, respectively.\n" %
        (len(data), len(data[0]), len(data[1]), len(data[-1])))
    return newdata

def removeEmptyRows(data):
    print("==> Removing empty rows from data... ", end="")
    cleandata=[]
    for m in range(len(data)):
        length = len(data[m])
        for n in range(length):
            if n < length-1 and data[m][n] != '':
                cleandata+=[data[m]]
                break
    print("done. \nData had %d empty rows, %d nonempty of %d total. \n" % 
        (len(data)-len(cleandata), len(cleandata), len(data)))
    return cleandata

def keepOnlyTheseColumns(data, colsToKeep):
    #colsToKeep should be the indicies of the desired columns
    #if left blank, keeps all cols
    if colsToKeep == []:
        print("==> Keeping all columns.\n")
        return data
    print("==> Keeping only these columns:")
    for n in range(len(colsToKeep)):
        print("Col "+str(colsToKeep[n])+": "+data[0][n])
    newdata=[]
    print("Working... ", end="")
    for n in range(len(data)):
        addline = []
        for col in colsToKeep:
            addline+=[data[n][col].replace('\n','')]
        newdata+=[addline]
    print("done. \n%d Columns remain of %d columns originally (in 1st line)" %
        (len(newdata[0]), len(data[0])))
    print("First line is now:")
    print(newdata[0])
    print("Second line is now:")
    print(newdata[1])
    print()
    return newdata

def rowLen(row):
    lastfullcol = -1
    for n in range(len(row)):
        if row[n]!='' and row[n]!='\n':
            lastfullcol = n
    return lastfullcol+1

def sanityCheck(data):
    #check 1: find the min and max length of the rows
    print("==> Sanity checking... ", end="")
    minlen=999999999999999999
    maxlen=0
    for n in range(len(data)):
        length = rowLen(data[n])
        if length > maxlen:
            maxlen = length
        if length < minlen:
            minlen = length
    print("done. \nThe number of columns per row ranges from %d to %d\n" %
        (minlen, maxlen))

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

def convertDataTypes(data, typelist):
    nSuccessfulTypes = 0
    nFailedTypes = 0
    if typelist==[]:
        return data
    #len(typelist) must be number of rows.

    print("==> Converting types. Upon error, leaving as string... ")
    newdata = []
    for n in range(len(data)):
        row = data[n]
        newrow = []
        if len(typelist) != len(data[n]):
            print("ERROR: length of typelist != length data[%d]." % n)
            sys.exit(1)

        for m in range(len(row)):
            typ = typelist[m].lower()
            typ = typ.replace(" ","")
            oldpt = row[m].lower()

            if typ in ["str", "string", "s"]:
                newrow+=[oldpt]
                nSuccessfulTypes+=1
                continue

            elif typ in ["float", "f"]:
                try:
                    newpt = float(oldpt)
                    newrow+=[newpt]
                    nSuccessfulTypes+=1
                    continue
                except ValueError:
                    nFailedTypes+=1
                    print("ConvertFail %d on row %d: cant cnvt '%s' to type '%s'." % 
                        (nFailedTypes, n, oldpt, typ))
                    newrow+=[oldpt]
                    continue

            elif typ in ["int", "i", "d"]:
                try:
                    newpt = int(oldpt)
                    newrow+=[newpt]
                    nSuccessfulTypes+=1
                    continue
                except ValueError:
                    nFailedTypes+=1
                    print("ConvertFail %d on row %d: cant cnvt '%s' to type '%s'." % 
                        (nFailedTypes, n, oldpt, typ))
                    newrow+=[oldpt]
                    continue

            elif typ in ["bool", "boolean", "b"]:
                truebools_lower = ["1", "true",  "t"]
                falsebools_lower= ["0", "false", "f"]
                if oldpt in truebools_lower:
                    newrow+=[True]
                    nSuccessfulTypes+=1
                    continue
                if oldpt in falsebools_lower:
                    newrow+=[False]
                    nSuccessfulTypes+=1
                    continue
                else:
                    nFailedTypes+=1
                    print("ConvertFail %d on row %d: cant cnvt '%s' to type '%s'." % 
                        (nFailedTypes, n, oldpt, typ))
                    newrow+=[oldpt]
                    continue

            elif typ in ["nparray", "np.array", "numpyarray", 
            "numpy.array", "np"]:
                try:
                    newpt = txtToNpArray(oldpt)
                    newrow+=[newpt]
                    nSuccessfulTypes+=1
                    continue
                except ValueError:
                    nFailedTypes+=1
                    print("ConvertFail %d on row %d: cant cnvt '%s' to type '%s'." % 
                        (nFailedTypes, n, oldpt, typ))
                    newrow+=[oldpt]
                    continue
            else:
                #this should never happen but just in case
                print("YOU SHOULD NEVER SEE THIS")
                newrow+=[oldpt]
                nFailedTypes+=1
        newdata+=[newrow]
    print("...done.\n")
    print("%d of %d (%d%%) type conversions were successful." % (nSuccessfulTypes,
        nSuccessfulTypes+nFailedTypes, 
        100*nSuccessfulTypes/(nSuccessfulTypes+nFailedTypes)))
    if nFailedTypes>0:
        print("Failed conversions were left as a string.")

    return newdata

def importAndClean(path, colsToKeep=[], typelist=[]):
    print("==============================="+"="*len(path))
    print("======IMPORTING/CLEANING %s======" % path)
    print("==============================="+"="*len(path))
    print("Importing... ",end="")
    newdata = []
    with open(path, "r") as f:
        reader = csv.reader(f, delimiter=",")
        for cell in reader:
            newdata += [cell]
        f.close()
    print("done.")
    print("First line is:\n", newdata[0])
    print("Second line is:\n", newdata[1])
    print()


    #check if inputted typelist is valid
    validtypes = ["string", "str", "s", 
    "float", "f", 
    "int", "i", "d", 
    "bool", "boolean","b", 
    "nparray", "np.array", "numpyarray", "numpy.array", "np"]

    typelist = [(x.lower()).replace(" ","") for x in typelist]
    for typ in typelist:
        if not (typ in validtypes):
            print("ERROR: You entered an invalid typename: '%s'." % typ)
            print("Valid type names are:")
            for n in range(len(validtypes)-1):
                print("%s, " % validtypes[n], end="")
            print(validtypes[-1], "\n")
            sys.exit(1)

    sanityCheck(newdata)
    newdata = keepOnlyTheseColumns(newdata, colsToKeep)
    newdata = convertDataTypes(newdata,typelist)
    #newdata = convertToNestedList(newdata) #deprecated: csvreader does it
    #newdata = removeEmptyRows(newdata)     #sometimes necessary
    return newdata

ans = input("""Ready to import from csv. 
    Have you backed up your csv file(s)? y/n: """).lower()
if not (ans in ["yes", "yy", "y"]):
    print("exiting. Please back up your data first!")
    sys.exit()
print()

#COLHEADER TYPES:
typelist=[
 "str",        "str",                            "np",         "f",               "f",         "f",                 "f",             "f",           "f",           "f",         "f",      "f",                      "s",                           "np",                "f",                      "f",                "np",             "f",                    "f",                  "f",                "f",                      "f",                        "f",                       "f",                    "f",             "f",                                 "f",                       "i",          "f",                     "np",                 "bool"]
#"ERROR CODE", "start time (YYYYMMDD-hhmmssms)", "v0 (kpc/s)", "mag(v0) (kpc/s)", "mag(v0/c)", "E0 (kg*kpc^2/s^2)", "theta0 (rads)", "phi0 (rads)", "q_b (A*kpc)", "mass (kg)", "dt (s)", "distance to halt (kpc)", "end time (YYYYMMDD-hhmmssms)","vel_final (kpc/s)", "mag(vel_final) (kpc/s)", "mag(vel_final/c)", "pos_final(kpc)", "mag(pos_final) (kpc)", "theta_final (rads)", "phi_final (rads)", "distance tracked (kpc)", "distance from start (kpc)","KE_final (kg*kpc^2/s^2)", "max velocity (kpc/s)", "max velocity/c","arclength in rho1kpc region (kpc)", "time after t=0 start (s)","iterations", "real runtime (real s)", "acc_final (kpc/s/s)","exit status (True/False)"]
#  rundata[0],                              [1],          [2],               [3],         [4],                 [5],             [6],           [7],           [8],         [9],     [10],                     [11],                           [12],               [13],                     [14],

filepath = "data/batchrun20160921-02395340.csv"
rundata = importAndClean(filepath, typelist=typelist)
colheaders = rundata[0]
rundata = rundata[1:]

# RUNDATA COLUMN HEADERS
# [0]  ERROR CODE
# [1]  start time (YYYYMMDD-hhmmssms)
# [2]  v0 (kpc/s)
# [3]  mag(v0) (kpc/s)
# [4]  mag(v0/c)
# [5]  E0 (kg*kpc^2/s^2)
# [6]  theta0 (rads)
# [7]  phi0 (rads)
# [8]  q_b (A*kpc)
# [9]  mass (kg)
# [10] dt (s)
# [11] distance to halt (kpc)
# [12] end time (YYYYMMDD-hhmmssms)
# [13] vel_final (kpc/s)
# [14] mag(vel_final) (kpc/s)
# [15] mag(vel_final/c)
# [16] pos_final(kpc)
# [17] mag(pos_final) (kpc)
# [18] theta_final (rads)
# [19] phi_final (rads)
# [20] distance tracked (kpc)
# [21] distance from start (kpc)
# [22] KE_final (kg*kpc^2/s^2)
# [23] max velocity (kpc/s)
# [24] max velocity/c
# [25] arclength in rho1kpc region (kpc)
# [26] time after t=0 start (s)
# [27] iterations
# [28] real runtime (real s)
# [29] acc_final (kpc/s/s)
# [30] exit status (True/False)

kgkpc2s2_To_GeV = 5.9427943e48

energydeltas = [transpose(rundata)[22][i] - transpose(rundata)[5][i] for 
    i in range(len(transpose(rundata)[22]))] #kg*kpc^2/s^2
energydeltas = [x*kgkpc2s2_To_GeV for x in energydeltas]

absenergydeltas = [abs(edelta) for edelta in energydeltas]

plt.hist(energydeltas)
# plt.hist(absenergydeltas)

plt.xlabel("Changes in monopole's kinetic energy, GeV")
plt.show()