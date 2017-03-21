from __future__ import print_function, division
try: input = raw_input
except NameError: pass

import csv
import sys
import numpy as np

# CSV import wrapper with type conversions and sanity checks and transpose function

# data = importcsv("ivs-dex.csv",
#     typelist=["int", "int", "int", "int"],
#     colsToKeep=[0,1,2,4])
# dataT = transpose(IVdata)

def transpose(data):
    trans = [list(x) for x in zip(*data)]
    return trans

def removeEmptyRows(data, print_to_console):
    # print("==> Removing empty rows from data... ", end="")
    cleandata=[]
    for m in range(len(data)):
        length = len(data[m])
        for n in range(length):
            if n < length-1 and data[m][n] != '':
                #if there is ANY nonempty stuff, add the whole row.
                cleandata+=[data[m]]
                break
    if print_to_console:
        print("Removed %d empty rows from data, of %d total. \n" % 
            (len(data)-len(cleandata), len(data)))
    return cleandata

def keepOnlyTheseColumns(data, colsToKeep, print_to_console):
    #colsToKeep should be the indicies of the desired columns
    #if left blank, keeps all cols
    if colsToKeep == []:
        if print_to_console:
            print("==> Keeping all columns.\n")
        return data
    if print_to_console:
        print("==> Keeping only these columns: ", end="")
    for n in range(len(colsToKeep)):
        print("Col " + str(colsToKeep[n]) + ": " + data[0][n] + ";  ", end="")
    newdata=[]
    for n in range(len(data)):
        addline = []
        for col in colsToKeep:
            addline+=[data[n][col].replace('\n','')]
        newdata+=[addline]
    if print_to_console:
        print("Done importing. First two lines are:\n", newdata[0])
        try: print(newdata[1])
        except IndexError: print()
        print()
    return newdata

def rowLen(row):
    lastfullcol = -1
    for n in range(len(row)):
        if row[n]!='' and row[n]!='\n':
            lastfullcol = n
    return lastfullcol+1

def sanityCheck(data, print_to_console):
    #check 1: find the min and max length of the rows
    minlen=999999999999999999
    maxlen=0
    for n in range(len(data)):
        length = rowLen(data[n])
        if length > maxlen:
            maxlen = length
        if length < minlen:
            minlen = length
    if print_to_console:
        if maxlen != minlen:
            print("Note: The number of columns per row ranges from %d to %d\n" % (minlen, maxlen))

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

def convertDataTypes(data, typelist, print_to_console, dont_convert_header):
    nSuccessfulTypes = 0
    nFailedTypes = 0
    if typelist==[]:
        return data
    #len(typelist) must be number of rows.
    if print_to_console:
        print("==> Converting types. Upon error, leaving as string... ")
    
    newdata = []
    rows_to_convert = range(len(data))
    if dont_convert_header: 
        newdata += [data[0]]
        rows_to_convert = range(1,len(data))
    
    for n in rows_to_convert:
        row = data[n]
        newrow = []
        if len(typelist) != len(data[n]):
            raise IndexError("length of typelist != length data[%d]." % n)

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
                    if n==0: continue
                    nFailedTypes+=1
                    if print_to_console:
                        print("ConvertFail %d on row %d col %d: cant cnvt '%s' to type '%s'." % 
                            (nFailedTypes, n, m, oldpt, typ))
                    newrow+=[oldpt]
                    continue

            elif typ in ["int", "i", "d"]:
                try:
                    newpt = int(oldpt)
                    newrow+=[newpt]
                    nSuccessfulTypes+=1
                    continue
                except ValueError:
                    if n==0: continue
                    nFailedTypes+=1
                    if print_to_console:
                        print("ConvertFail %d on row %d col %d: cant cnvt '%s' to type '%s'." % 
                            (nFailedTypes, n, m, oldpt, typ))
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
                    if n==0: continue
                    nFailedTypes+=1
                    if print_to_console:
                        print("ConvertFail %d on row %d col %d: cant cnvt '%s' to type '%s'." % 
                            (nFailedTypes, n, m, oldpt, typ))
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
                    if n==0: continue
                    nFailedTypes+=1
                    if print_to_console:
                        print("ConvertFail %d on row %d col %d: cant cnvt '%s' to type '%s'." % 
                            (nFailedTypes, n, m, oldpt, typ))
                    newrow+=[oldpt]
                    continue
            else:
                #this should never happen but just in case
                print("YOU SHOULD NEVER SEE THIS")
                newrow+=[oldpt]
                nFailedTypes+=1
        newdata+=[newrow]

    if print_to_console:    
        print("...done.\n")
        print("%d of %d (%.3f%%) type conversions were unsuccessful. " % (nFailedTypes,
            nSuccessfulTypes+nFailedTypes, 
            100*nFailedTypes/(nSuccessfulTypes+nFailedTypes)), end="")
        if nFailedTypes>0:
            print("Failed conversions were left as a string.")

    return newdata

def importcsv(path, colsToKeep=[], typelist=[], print_to_console=True, dont_convert_header=True,
    remove_empty_rows=True):
    if print_to_console:
        ans = input("Ready to import from csv. Have you backed up your csv file(s)? y/n: ")
        if not (ans in ["yes", "yy", "y"]):
            print("exiting. Please back up your data first!")
            sys.exit()
        print()

    newdata = []
    with open(path, "r") as f:
        reader = csv.reader(f, delimiter=",")
        for cell in reader:
            newdata += [cell]
        f.close()
    if print_to_console:
        print("Done importing. First two lines are:\n", newdata[0])
        try: print(newdata[1])
        except IndexError: print()
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
            raise NameError( ("You entered an invalid typename in typelist '%s'. \n" % typ) + 
                "Please redo typelist.")

    sanityCheck(newdata, print_to_console)
    newdata = keepOnlyTheseColumns(newdata, colsToKeep, print_to_console)
    newdata = convertDataTypes(newdata, typelist, print_to_console, dont_convert_header)
    if remove_empty_rows:
        newdata = removeEmptyRows(newdata, print_to_console)     #sometimes necessary
    return newdata

# data = importcsv("ivs-dex.csv",
#     typelist=["int", "int", "int", "int"],
#     colsToKeep=[0,1,2,4])
# dataT = transpose(IVdata)