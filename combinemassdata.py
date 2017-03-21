from __future__ import division, absolute_import, print_function, unicode_literals
try: input = raw_input
except NameError: pass
try: range = xrange
except NameError: pass

import csvimport as cv
import sys
import os
import csv

# EXTRACTING THE DATA FROM FILES
filelist = os.listdir("data")
answer = raw_input("Is the data folder backed up? y/n: ")
if answer == "no" or answer == "n" or answer == "N":
    print("exiting.")
    sys.exit()
print("Good to know. Combining file number 1 of %d" % len(filelist))

data_from_csvs = []
first_import = cv.importcsv("data/" + str(filelist[0]), print_to_console = False, remove_empty_rows=False)

if first_import[0][0] != "SETTINGS:": settings_rownum = first_import.index(["SETTINGS:"])
else: settings_rownum = 0

if first_import[4] != "RUN DATA:": rundata_rownum = first_import.index(["RUN DATA:"])
else: rundata_rownum = 4

# these settings should be the same across all files
st = first_import[settings_rownum + 2]
timestamp, dt, q_b = first_import[rundata_rownum+2][2:5]

runparameters = [dt, q_b]

# do some sanity checks
st = map(lambda s: str(s), st)
if (st[0] != "True" or st[7] != "1" or st[9] != "True" or st[10] != "True" or st[11] != "True"):
    print("SETTINGS WERE BAD FOR FIRST CSV! EXITING.\n")
    sys.exit(1)

combined_data = []
timestamp_list = []

for n in range(1, len(filelist)):
    print("\rGood to know. Combining file number %d of %d" % (n+1, len(filelist)))
    data = cv.importcsv("data/" + str(filelist[n]), print_to_console = False, remove_empty_rows=False)
    st_ = data[settings_rownum+2]
    if st_ != st:
        print("SETTINGS WERE DIFFERENT IN FILE NUMBER %d (numbered starting from 0)" % n)
        sys.exit(1)

    dt_, q_b_ = data[rundata_rownum+2][3:5]
    timestamp_list += data[rundata_rownum+2][2]

    if [dt_, q_b_] != runparameters:
        print("PARAMETERS WERE DIFFERENT IN FILE NUMBER %d (numbered starting from 0)" % n)
        print("First file parameters dt, q_b were", runparameters)
        print("This file parameters are", [dt_, q_b_])
        sys.exit(1)

    combined_data += data[rundata_rownum+2:]

timestamp_list = sorted(timestamp_list)

# add the settings to combined data
output_csv_header = []
output_csv_header += [["COMBINED DATA FOR %d DIFFERENT MASSES" % len(filelist)]]
output_csv_header += [["TIMESTAMPS FROM, TO:", "",timestamp_list[0], timestamp_list[-1]]]
output_csv_header += [[]]
output_csv_header += first_import[settings_rownum:settings_rownum+3]
output_csv_header += [[]]
output_csv_header += first_import[rundata_rownum:rundata_rownum+2]

output_csv = output_csv_header + combined_data

outputfilename = "AGGREGATEDbatchrun" + timestamp + ".csv"
print("Writing to file %s ... " % outputfilename, end="")
with open(outputfilename, 'wb') as file:
    mywriter = csv.writer(file)
    for line in output_csv:
        mywriter.writerow(line)
    file.close()
print("done.\n")