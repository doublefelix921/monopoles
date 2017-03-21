#take 5 args and print them 10x each, running in parallel

import sys
import os

from time import strftime as tstr

os.chdir(os.path.dirname(sys.argv[0]))

blah = open("testrun"+tstr("%Y.%m.%d.%H.%M.%S")+".txt","a")

blah.write(str(sys.argv[1])+"\n")
blah.write(str(sys.argv[2])+"\n")
blah.write(str(sys.argv[3])+"\n")
blah.write(str(sys.argv[4])+"\n")
blah.write(str(sys.argv[5])+"\n")


blah.close()

home/feist/riatest